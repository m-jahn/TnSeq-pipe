#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for NewerThan()

my $nCPUs = $ENV{MC_CORES} || (`egrep -c '^processor' /proc/cpuinfo`);
my $maxload = $ENV{SUBMITTER_MAXLOAD} || $ENV{MC_CORES} || $nCPUs+1;
my $nUse = $ENV{SUBMITTER_NJOBS} || $ENV{MC_CORES} || int(0.6 * ($nCPUs+1));
$nUse = 1 if $nUse < 1;

my $usage = "Usage: RunRapSearch.pl [ -test ] [ -nCPU $nUse ] nickname1 ... nicknameN";
my $test;
die $usage unless GetOptions('test' => \$test, 'nCPU=i' => \$nUse);
my $rapsearch = "$Bin/rapsearch";
die "No such executable: $rapsearch\n" unless -x $rapsearch;
my @dbs = ("sprot","KEGG","metacyc");
foreach my $db (@dbs) {
  my $rap = lc($db) . ".rap";
  die "No such file: $rap\n" unless -e $rap;
}
my $dir = "g/blast_results";
die "No such directory: $dir\n" unless -d $dir;

my @orgs = @ARGV;
die $usage if @orgs == 0;

my @work = (); # list of organism to rebuild, and either sprot or KEGG
foreach my $org (@orgs) {
    die "No such directory: g/$org" unless -d "g/$org";
    my $seqfile = "g/$org/aaseq2";
    die "No such file: $seqfile (you may need to run MakeBBH.pl or SetupDoms.pl)"
        unless -e $seqfile;
    foreach my $db (@dbs) {
        my $rapfile = lc($db) . ".rap";
        my $hits = "$dir/${db}_${org}.m8";
        unless (NewerThan($hits, $seqfile)
                && NewerThan($hits, $rapfile)) {
            print STDERR "Rerunning $org vs. $db\n";
            push @work, [ $org, $db ];
        }
    }
}

print STDERR "Running rapsearch with $nUse cpus\n" if @work > 0;
if (defined $test) {
    print STDERR "Testing mode, " . scalar(@work) . " jobs will not actually be run\n";
} else {
    foreach my $row (@work) {
        my ($org,$db) = @$row;
        my $rapfile = lc($db) . ".rap";
        my $outpre = "$dir/${db}_${org}";

        system("feba/bin/rapsearch -q g/$org/aaseq2 -t a -e -2 -b 30 -v 30 -z $nUse -d $rapfile -o $outpre") == 0
            || die "rapsearch failed to for $org with database $rapfile";
    }
}
