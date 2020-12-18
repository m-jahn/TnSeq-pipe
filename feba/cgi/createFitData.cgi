#!/usr/bin/perl -w

#######################################################
## createExpData.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page generates output to STDOUT and presents it as
# a download, which includes all of the experimental fitness data
# for a single organism. It is run from org.cgi. By default
# it makes a file of fitness values, but it can
# save t values instead.
#
#
# Required CGI parameters:
# orgId -- which organism to search for
#
# Optional CGI parameters:
# expName -- which experiments to include. If not specified, include all experiments.
# t -- report t values rather than fitness values

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi = CGI->new;

my $orgId = $cgi->param('orgId');
die "Must specify orgId" unless defined $orgId;
my @expNames = $cgi->param('expName');
my $show_t = $cgi->param('t');

my $dbh = Utils::get_dbh();

# gather all of the data needed

my $gene = $dbh->selectall_arrayref("SELECT orgId, locusId, sysName, gene, desc FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

my $exp = $dbh->selectall_arrayref("SELECT expName, expDesc FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);
my %exp = map { $_->{expName} => $_ } @$exp;
foreach my $expName (@expNames) {
  die "Invalid expName $expName in org $orgId\n"
    unless exists $exp{$expName};
}

my $fetchname = $show_t ? "t" : "fit";
my $fitness = $dbh->selectall_arrayref("SELECT locusId, expName, $fetchname FROM GeneFitness WHERE orgId = ?", {}, $orgId);

# This may store t values rather than actual fitness values
my %fit=();
foreach my $frow(@$fitness) {
	my ($locusId, $expName, $fit) = @$frow;
	$fit{$locusId}{$expName} = $fit;
}

my $prefix = $show_t ? "t" : "fit";
my @pieces = ($prefix, "organism", $orgId);
push @pieces, "selected" if @expNames > 0;
my $outfile = join("_", @pieces) . ".tsv";
print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=$outfile\n\n";

@expNames = map { $_->{expName} } @$exp
  if @expNames == 0;
my @expsUse = map $exp{$_}, @expNames;

# print the header row
print join("\t",
           qw{orgId locusId sysName geneName desc},
           map { $_->{expName} . " " . $_->{expDesc} } @expsUse) . "\n";

# print the data row by row
foreach my $row (@$gene) {
  my $locusId = $row->{locusId};
  next if !exists $fit{$locusId};
  my @out = ($row->{orgId}, $locusId, $row->{sysName}, $row->{gene}, $row->{desc});
  push @out, map { sprintf("%.3f", $fit{$locusId}{$_}) } @expNames;
  print join("\t", @out) . "\n";
}
