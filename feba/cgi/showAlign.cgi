#!/usr/bin/perl -w
#######################################################
## showAlign.cgi
##
## Copyright (c) 2017 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# query (orgId:locusId or the defline)
# subject (orgId:locusId)
#
# Optional CGI parameters:
# querySequence (if query is not from the database)

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use Bio::SeqIO;

use lib "../lib";
use Utils;
sub ProtInfo($); # orginfo:locusId => HTML text

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information
my $query = $cgi->param('query') || "";
my $subject = $cgi->param('subject') || "";
my $querySequence = $cgi->param('querySequence');
die "Must specify query" unless $query;
die "Must specify subject" unless $subject;
die "Invalid subject" unless $subject =~ m/:/;
die "Invalid query" unless defined $querySequence || $query =~ m/:/;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $myDB = Utils::blast_db();
my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $filename = $procId . $timestamp;
my $tmpDir = Utils::tmp_dir();
my $prefix = "$tmpDir/$filename";
my $file1 = "$prefix.1";
my $file2 = "$prefix.2";
my $alnfile = "$prefix.aln";

my $fastacmd = '../bin/blast/fastacmd';
die "No such command: $fastacmd" unless -x $fastacmd;

my $bl2cmd = '../bin/blast/bl2seq';
die "No such command: $bl2cmd" unless -x $bl2cmd;

print
  header(),
  Utils::start_page("Pairwise Alignments"),
  h2("Pairwise Alignments");

# set up the first sequence
if (defined $querySequence) {
  open(FILE, ">", $file1) || die "Cannot write to $file1";
  # do not use actual defline, could cause bl2seq to fail
  print FILE ">query\n$querySequence\n";
  close(FILE) || die "Error writing to $file1";
} else {
  system($fastacmd, "-d", $myDB, "-s", $query, "-o", $file1) == 0
    || die "$fastacmd for $query failed: $!";
}
system($fastacmd, "-d", $myDB, "-s", $subject, "-o", $file2) == 0
    || die "$fastacmd for $query failed: $!";

system($bl2cmd,
       "-p", "blastp",
       "-e", 0.01,
       "-F", "m S",
       "-i", $file1,
       "-j", $file2,
       "-o", $alnfile) == 0
  || die "$bl2cmd failed on $file1 $file2";
my @lines = ();
open(IN, "<", $alnfile) || die "Cannot read $alnfile";
@lines = <IN>;
@lines = map { chomp; $_; } @lines;
close(IN) || die "Error reading $alnfile";

my $qlen = 0;
if (defined $querySequence) {
  $qlen = length($querySequence);
} else {
  $qlen = length( Bio::SeqIO->new(-file => $file1, -format => 'fasta')->next_seq()->seq );
}
my $slen = length( Bio::SeqIO->new(-file => $file2, -format => 'fasta')->next_seq()->seq );

unlink($file1);
unlink($file2);
unlink($alnfile);

# remove lines up to the first alignment ("Score")
do {
    shift @lines
} until $lines[0] =~ m/Score/i || $lines[0] =~ m/No hits/i || @lines == 1;
my @out = ();
foreach my $line (@lines) {
    if ($line =~ m/^Lambda\s+K\s+H$/) {
        last;
    } else {
        push @out, $line;
    }
}

print
  p("Query, $qlen a.a.,",
    defined $querySequence ? $query : &ProtInfo($query)),
  p("Subject, $slen a.a.,",
    &ProtInfo($subject)),
  pre(join("\n", @out));
Utils::endHtml($cgi);
exit(0);

sub ProtInfo($) {
  my ($spec) = @_;
  my ($orgId, $locusId) = split /:/, $spec;
  die "Invalid spec $spec" unless $orgId && $locusId;
  my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
  die "Unknown gene $spec" unless defined $gene;
  my $showName = $gene->{sysName} || $gene->{locusId};
  $showName .= " " . $gene->{gene} if $gene->{gene};
  return join(" ",
              Utils::gene_link($dbh, $gene, "desc", "domains.cgi"),
              "from",
              a( { href => "org.cgi?orgId=$orgId" },
                 $orginfo->{$orgId}{genome} ));
}
