#!/usr/bin/perl -w
#######################################################
## getNtSeq.cgi - fetch nucleotide sequence of a
##                scaffold or around a gene
##
## Copyright (c) 2018 University of California
##
## Author: Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
#
# Optional CGI parameters:
# locusId -- which gene to fetch sequence around
# scaffoldId -- which scaffold to fetch the sequence of
# begin,end -- which part of the scaffold to fetch
#	end < begin means use - strand

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
use Utils;
use FEBA_Utils; # for reverseComplement
sub PrintFasta($$);

my $cgi=CGI->new;
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";

print header(-type => 'text/plain');

if (my $locusId = $cgi->param('locusId')) {
  my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
                                   {}, $orgId, $locusId);
  die "No such gene $locusId in $orgId" unless $gene;
  my $scaffoldId = $gene->{scaffoldId};
  my $seq = $dbh->selectrow_array("SELECT sequence FROM ScaffoldSeq WHERE orgId = ? AND scaffoldId = ?",
                                  {}, $orgId, $gene->{scaffoldId});
  my $begin = $gene->{begin};
  my $end = $gene->{end};
  die "begin must be less than end" unless $begin < $end;
  die "Invalid begin/end" unless $begin >= 1 && $end <= length($seq);
  # Logic as if gene is on + strand
  my $geneseq = substr($seq, $begin-1, $end-$begin+1);
  my $beginUp = $begin - 200;
  $beginUp = 1 if $beginUp < 1;
  my $upseq = substr($seq, $beginUp-1, $begin-$beginUp);
  my $endDn = $end + 200;
  $endDn = length($seq) if $endDn > length($seq);
  my $dnseq = substr($seq, $end, $endDn - $end);
  # Transform if on - strand
  if ($gene->{strand} eq "-") {
    ($upseq,$geneseq,$dnseq) = (reverseComplement($dnseq), reverseComplement($geneseq), reverseComplement($upseq));
  }
  my $showGene = $gene->{sysName} || $gene->{locusId};
  my $uplen = length($upseq);
  my $dnlen = length($dnseq);
  PrintFasta("$showGene:upstream$uplen from $genome", $upseq);
  PrintFasta("$showGene:coding from $genome", $geneseq);
  PrintFasta("$showGene:downstream$dnlen from $genome", $dnseq);
} elsif (my $scaffoldId = $cgi->param('scaffoldId')) {
  my $seq = $dbh->selectrow_array("SELECT sequence FROM ScaffoldSeq WHERE orgId = ? AND scaffoldId = ?",
                                  {}, $orgId, $scaffoldId);
  die "No sequence for $scaffoldId in $genome\n" unless $seq;
  my $begin = $cgi->param('begin');
  my $end = $cgi->param('end');
  if ($begin && $end) {
    $begin = 1 if $begin < 1;
    $end = 1 if $end < 1;
    $begin = length($seq) if $begin > length($seq);
    $end = length($seq) if $end > length($seq);
    my ($strand, $subseq);
      if ($begin <= $end) {
        $strand = "+";
        $subseq = substr($seq, $begin-1, $end-$begin+1);
      } else {
        $strand = "-";
        $subseq = reverseComplement(substr($seq, $end-1, $begin-$end+1));
      }
    PrintFasta("$scaffoldId:$begin:$end:$strand $orgId $genome", $subseq);
  } else {
    PrintFasta("$scaffoldId $orgId $genome", $seq);
  }
} else {
  print "No locus or scaffold specified\n";
}

sub PrintFasta($$) {
  my ($header, $seq) = @_;
  print ">$header\n";
  my $sz = 60;
  for (my $i = 0; $i < length($seq); $i += $sz) {
    print substr($seq, $i, $sz)."\n";
  }
}
