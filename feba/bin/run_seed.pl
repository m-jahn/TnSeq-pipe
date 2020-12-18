#!/usr/bin/perl -w
# Use the Network-based SEED API to annotation the proteins in the specified nickname
# Assumes that the seed API software
#	http://blog.theseed.org/downloads/sas.tgz
# is installed in sas/
# This script is modified from server_paper_example6.pl
#	http://servers.nmpdr.org/sapling/server.cgi?code=server_paper_example6.pl
use strict;
use lib "sas/lib";
use lib "sas/modules/lib";

use ANNOserver;
use Data::Dumper;

die "Usage: run_seed.pl nickname\n" unless @ARGV==1;
my $dir = $ARGV[0];
die "Invalid nickname: $dir\n" unless -d "g/$dir";
my $infile = "g/$dir/aaseq";
die "No aaseq file in g/$dir/\n" unless -e $infile;
my $annoObject = ANNOserver->new();
my $blastOption = 0;
open(IN, "<", $infile) || die "Cannot read $infile";
my $outfile = "g/$dir/seedanno.tab";
my $tmpout = "$outfile.tmp.$$";
open(OUT, ">", $tmpout) || die "Cannot write to $tmpout";
print STDERR "Running SEED server for nickname $dir\n";
my $results = $annoObject->assign_function_to_prot(-input => \*IN,
                                                   -kmer => 8,
                                                   -scoreThreshold => 3,
                                                   -assignToAll => $blastOption);
my $n = 0;
my $nTot = 0;
while (my $data = $results->get_next()) {
    my ($id, $role, $genomeSet, undef, $hits) = @$data;
    $nTot++;
    # Only proceed if a role was found.
    if ($role) {
        print OUT join("\t", $id, $role) . "\n";
        $n++;
    }
}
close(IN) || die "Error reading $infile";
close(OUT) || die "Error writing to $tmpout";
rename($tmpout,$outfile) || die "Error renaming $tmpout to $outfile";
print STDERR "Wrote SEED annotations for $n of $nTot to $outfile\n";

