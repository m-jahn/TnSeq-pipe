#!/usr/bin/perl -w

#######################################################
## orgSeqs.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# This page generates an amino acid fasta sequence file
# or nucleotide (genome) fasta file for an organism.
#
# Required CGI parameters:
# orgId -- which organism to search for
# optional -- type=nt for genome sequence

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi = CGI->new;

my $orgId = $cgi->param('orgId');
die "No orgId" if !defined $orgId || $orgId eq "";

my $type = $cgi->param('type') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};

my $suffix = $type eq "nt" ? "fna" : "faa";

print "Content-Type:application/x-download\n";
print "Content-Disposition: attachment; filename=organism_$orgId.$suffix\n\n";

if ($type eq "nt") {
  my $sc = $dbh->selectall_arrayref("SELECT * from ScaffoldSeq WHERE orgId = ?",
                                   { Slice => {} }, $orgId);
  foreach my $row (@$sc) {
    print ">" . $row->{scaffoldId} . "\n";
    my @whole = $row->{sequence} =~ /.{1,60}/g;
    map { print "$_\n" } @whole;
  }
} else {
  my $genes = $dbh->selectall_hashref("SELECT * FROM Gene WHERE orgId = ?",
                                      "locusId", # field to index hash by
                                      {Slice => {}}, $orgId);
  die "No genes in $orgId" if scalar(keys %$genes) == 0;
  my $aaseqs = Utils::blast_db();
  open(FAA, "<", $aaseqs) || die "Cannot read $aaseqs";
  my $printing = 0;
  while(<FAA>) {
    chomp;
    if (m/^>(.*)$/) {
      my ($orgLocus,$locusId) = split /:/, $1;
      if ($orgLocus eq $orgId) {
        $printing = 1;
        my $gene = $genes->{$locusId};
        die "Unrecognized gene $locusId in $orgId" if !defined $gene;
        print ">$orgId:$locusId $gene->{sysName} $gene->{desc}\n";
      } else {
          $printing = 0;
        }
    } else {
      print $_."\n" if $printing;
    }
  }
  close(FAA) || die "Error reading $aaseqs";
}
