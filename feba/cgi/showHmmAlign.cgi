#!/usr/bin/perl -w
#######################################################
## showAlign.cgi
##
## Copyright (c) 2020 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId, locusId,
# domainDb ('PFam' or 'TIGRFam'),
# and domainName (such as 'TIGR01926' or 'Cupin_6')

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use Bio::SeqIO;
use IO::Handle; # for autoflush

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information
my $orgId = $cgi->param('orgId');
my $locusId = $cgi->param('locusId');
my $domainDb = $cgi->param('domainDb');
my $domainName = $cgi->param('domainName');
die "Must specify orgId, locusId, domainDb, domainName"
    unless $orgId && $locusId && $domainDb && $domainName;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $myDB = Utils::blast_db();
my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $filename = $procId . $timestamp;
my $tmpDir = Utils::tmp_dir();
my $prefix = "$tmpDir/$filename";

die "Unknown orgId" unless exists $orginfo->{$orgId};
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
die "Unknown locusId" unless defined $gene->{locusId};

die "Invalid domainDb" unless $domainName =~ m/^[a-zA-Z0-9._-]+$/;
die "Invalid domainName" unless $domainName =~ m/^[a-zA-Z0-9._-]+$/;

my $fastacmd = '../bin/blast/fastacmd';
die "No such command: $fastacmd" unless -x $fastacmd;

my $hmmFile = "$prefix.hmm";
my $hmmDb = "../hmm/" . lc($domainDb) . ".hmm";
die "No such file: $hmmDb" unless -e $hmmDb;
system("../bin/hmmfetch $hmmDb $domainName > $hmmFile") == 0
  || die "hmmfetch failed -- $domainName not a valid domain name in $hmmDb";
my $hmmLen;
open(my $fhHmm, "<", $hmmFile) || die "Cannot read $hmmFile";
while (my $line = <$fhHmm>) {
  chomp $line;
  $hmmLen = $1 if $line =~ m/^LENG[\t ]+(\d+)$/;
}
close($fhHmm) || die "Error reading $hmmFile";

my $faaFile = "$prefix.faa";
system($fastacmd, "-d", $myDB, "-s", "${orgId}:${locusId}", "-o", $faaFile) == 0
    || die "$fastacmd for ${orgId}:${locusId} failed: $!";
my $in = Bio::SeqIO->new(-file => $faaFile,-format => 'fasta');
my $seq = $in->next_seq()->seq;

print
  header(),
  Utils::start_page("HMMer Alignments for $domainName"),
  h2("HMMer Alignments");

print p("Aligning the $domainName family to the amino acid sequence of",
        Utils::gene_link($dbh, $gene, 'desc', 'domains.cgi')),
  p("HMM length:", $hmmLen, "amino acids"),
  p("Sequence length:", length($seq), "amino acids");


autoflush STDOUT 1; # so preliminary results appear
print "\n";

print "<pre>\n";
system("../bin/hmmsearch", $hmmFile, $faaFile) == 0
  || die "hmmsearch failed: $!";
print "</pre>\n";
unlink($faaFile);
unlink($hmmFile);

print p(small(qq{The lines below each alignment block (ending in "PP") shows the
                 posterior probability that the positions align.
                 This is roughly the expected accuracy.
                 "*" means at least 95% accuracy,
                 0 means at most 5% accuracy, and 1-9 means about 10-90% confidence.}));

my %SSNames = ("B" => "beta-bridge (isolated)",
               "C" => "coil/loop",
               "E" => "beta-strand (hydrogen-bonded extended strand)",
               "G" => "3/10-helix",
               "H" => "alpha-helix",
               "I" => "pi (5)-helix",
               "S" => "bend",
               "T" => "turn (hydrogen-bonded)");

print p(small(qq{If the consensus structures are noted above the HMM sequence (lines ending in "CS"),
                 then the codes are:}));
print "<UL>\n";
foreach my $ss (sort keys %SSNames) {
  print li(small("$ss: $SSNames{$ss}"));
}
print "</UL>\n";

Utils::endHtml($cgi);
exit(0);
