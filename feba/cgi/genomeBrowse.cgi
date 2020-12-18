#!/usr/bin/perl -w
use strict;
#######################################################
## genomeBrowse.cgi -- show genes in a region of a
## scaffold, along with additional user-defined objects
##
## Copyright (c) 2018 University of California
##
## Author: Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# scaffoldId -- which scaffold
#
# Optional CGI arguments:
# object -- 1 or more additional objects, specified by a colon-delimited set of arguments such as
#	b:1000:e:2000:s:1:n:name:d:popup
#	begin, end, and name must be specified, and begin should be less than end
#	object is shown as unstranded if strand is absent
#	use 1 for + strand, not "+", because of CGI encoding issues
#	d is shown as popup text and is optional
# begin, end -- which part of the scaffold to show (optional)
# zoom = "in" or "out"
# pan = "left" or "right"
#
# begin and end must be specified unless objects are given

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $dbh = Utils::get_dbh();

my $orgId = $cgi->param('orgId') || die "No orgId found";
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Unknown orgId $orgId\n";
my $scaffoldId = $cgi->param('scaffoldId') || die "Invalid scaffoldId";
my $begin = $cgi->param('begin');
my $end = $cgi->param('end');

my %keynames = ("b" => "begin", "e" => "end", "s" => "strand", "n" => "name", "d" => "desc");
my @objects = ();
foreach my $objspec ($cgi->param('object')) {
  my @parts = split /:/, $objspec;
  my $row = { 'object' => 1 };
  while (@parts > 0) {
    my $key = shift @parts;
    die "Unknown key $key\n" unless exists $keynames{$key};
    my $value = shift @parts;
    die "Wrong number of fields in $objspec\n" unless defined $value;
    $row->{$keynames{$key}} = $value;
  }
  die "Invalid object $objspec\n"
    unless exists $row->{begin} && exists $row->{end} && exists $row->{name};
  # To build the URL for getNtSeq.cgi, may need to reverse begin and end
  my ($begin, $end) = ($row->{begin}, $row->{end});
  ($begin,$end) = ($end,$begin)
    if exists $row->{strand} && ($row->{strand} eq "-" || $row->{strand} eq "-1");
  $row->{URL} = "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end";
  push @objects, $row;
}
@objects = sort { $a->{begin} <=> $b->{begin} } @objects;

if (@objects > 0 && (!defined $begin || !defined $end)) {
  $begin = $objects[0]{begin} - 100;
  $end = $objects[-1]{end} + 100;
  if ($begin - $end + 1 < 4000) {
    my $center = int(($begin+$end)/2);
    $begin = $center - 2000;
    $end = $center + 2000;
  }
}

my ($sclen) = $dbh->selectrow_array("SELECT length(sequence) FROM ScaffoldSeq WHERE orgId = ? AND scaffoldId = ?",
                                    {}, $orgId, $scaffoldId);
die "Invalid scaffold $scaffoldId" unless $sclen;

$begin = 1 if $begin < 1;
if (!defined $end) {
  $end = $begin + 8000;
  $end = $sclen if $end > $sclen;
  $begin = $sclen - 8000 if $begin > $end-4000;
  $begin = 1 if $begin < 1;
}
die "End is too low\n" if $end < 1;
$begin < $end || die "begin must be less than end\n";

unless ($end <= $sclen) {
  print header;
  Utils::fail($cgi, "End is past the end of the scaffold");
}

my $zoom = $cgi->param('zoom') || "";
my $initwidth = $end - $begin + 1;
if ($zoom eq "in") {
    $begin += 0.2 * $initwidth;
    $end -= 0.2 * $initwidth;
} elsif ($zoom eq "out") {
    $begin -= 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
my $pan = $cgi->param('pan') || "";
if ($pan eq "left") {
    $begin -= 0.4 * $initwidth;
    $end -= 0.4 * $initwidth;
} elsif ($pan eq "right") {
    $begin += 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
$begin = int($begin);
$end = int($end);
# prevent panning off the edge
$begin = 1 if $begin < 1;
$begin = $sclen if $begin > $sclen;
$end = 1 if $end < 1;
$end = $sclen if $end > $sclen;

if ($end-$begin+1 < 100) {
  print header;
  Utils::fail($cgi, "Too small a region to show");
}
if ($end-$begin+1 > 50*1000) {
  print header;
  Utils::fail($cgi, "Too large a region to show");
}
my $title = "Genes";
my $nobj = scalar(@objects);
$title .= " (and $nobj objects)" if $nobj;

my $scaffoldShow = $scaffoldId;
print header,
  Utils::start_page("Browse $genome"),
  q{<div id="ntcontent">},
  h2("Browse", a({-href => "org.cgi?orgId=$orgId"}, $genome)),
  p($title, "for nucleotides", Utils::commify($begin), "to",
    Utils::commify($end), "on scaffold $scaffoldId:");

my $genes = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE orgId = ? AND scaffoldId = ? AND Gene.end >= ? AND Gene.begin <= ? ORDER by Gene.begin",
                                     { Slice => {} }, $orgId, $scaffoldId, $begin, $end);
my @show = @$genes;
push @show, @objects;
@show = sort { $a->{begin} <=> $b->{begin} } @show;

# specify no locus to highlight
print Utils::geneArrows(\@show, undef, $begin, $end);

my @objhidden = ();
foreach my $objspec ($cgi->param('object')) {
  push @objhidden, hidden(-name => 'object', -value => $objspec, -override => 1);
}

print
  start_form(-name => 'input', -method => 'GET', -action => 'genomeBrowse.cgi'),
  hidden( -name => 'orgId', -value => $orgId, -override => 1),
  hidden( -name => 'scaffoldId', -value => $scaffoldId, -override => 1),
  hidden( -name => 'begin', -value => $begin, -override => 1),
  hidden( -name => 'end', -value => $end, -override => 1),
  join("\n", @objhidden),
  p({-class => "buttons", style=>"max-width:500px; line-height:40px; white-space:nowrap;"},
    "Zoom:", submit('zoom','in'), submit('zoom','out'), "\tPan:", submit('pan','left'), submit('pan','right')),
  end_form;

print
  start_form(-name => 'input', -method => 'GET', -action => 'strainTable.cgi'),
  hidden( -name => 'orgId', -value => $orgId, -override => 1),
  hidden( -name => 'scaffoldId', -value => $scaffoldId, -override => 1),
  hidden( -name => 'begin', -value => $begin, -override => 1),
  hidden( -name => 'end', -value => $end, -override => 1),
  join("\n", @objhidden),
  p({-class => "buttons", style=>"margin-top: 2em; text-align: left; white-space:nowrap; line-height:40px;"}, "Add experiment(s): ",
    textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
    submit('Add','Add')),
  end_form;

print br(),
  p("Or see this region's",
    a({ -href => "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end" },
      "nucleotide sequence"),
    "or the entire sequence of scaffold",
    a({ -href => "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId" },
      $scaffoldId),
    "or the",
    a({ -href => "https://iseq.lbl.gov/genomes/genome/find?name=$orgId" }, "ENIGMA genome browser"));

print q{</div>};
$dbh->disconnect();
Utils::endHtml($cgi);
