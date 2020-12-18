#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for NewerThan
sub AddOrgName($$$); # Convert aaseq file to aaseq2

{
    my $blastdir = "$Bin/blast";
    my $outdir = ".";
    my $datadir = "html"; # where to find the fit_logratios_good.tab files
    my $gdir = "g"; # where to find the amino acid sequence files
    my $nCPU = undef;
    my $usage = <<END
Usage: MakeBBH.pl [ -datadir $datadir]
                  [ -gdir $gdir ]
		  [ -nCPU nCPU]
                  [ -outdir $outdir ]
                  nickname1 nickname2 ... nicknameN

MakeBBH.pl runs BLASTp and then computes bidirectional best hits and
paralogs (using bbh.pl and paralogs.pl). It uses identifiers of the
form organismId:locusId.

Requirements:

For each organism, g/nickname/aaseq should have the protein sequences
for that organism (in fasta format). aaseq2 in that directory should
have the sequences with organismid:locusId identifiers; if not,
MakeBBH.pl will create that file.

If html/nickname/fit_logratios_good.tab does not exist, it issues a
warning.

The BLAST executables (formatdb and blastall) should be in the
feba/bin/blast/ directory.

To avoid recomputing all pairs when adding another genome,
stores the results for each pair of genomes in
g/blast_results/blastp_nickname1_vs_nickname2

Results:

outdir/aaseqs.bbh -- table of all orthologs for each gene, with 1 row per gene
outdir/aaseqs.scores -- one row per pair of orthologs
outdir/aaseqs.para -- one row per pair of paralogs

(Also uses outdir/aaseqs.blastp as a temporary file)
END
    ;

    die $usage unless GetOptions('datadir=s' => \$datadir,
                                 'outdir=s' => \$outdir,
				 'gdir=s' => \$gdir,
				 'nCPU=i' => \$nCPU );
    my @orgs = @ARGV;
    die "MakeBBH.pl must have two or more organisms\n" unless scalar(@orgs) >= 2;
    die "No such directory: $blastdir" unless -d $blastdir;
    die "No such executable: $blastdir/formatdb" unless -x "$blastdir/formatdb";
    die "No such executable: $blastdir/blastall" unless -x "$blastdir/blastall";

    # Check for inputs and make the aaseq2 files
    my @seq2files = ();
    foreach my $org (@orgs) {
	my $seqfile = "$gdir/$org/aaseq";
	die "No aaseq file for $org: $seqfile" unless -e $seqfile;
	my $fitfile = "$datadir/$org/fit_logratios_good.tab";
	if (! -e $fitfile) {
	    print STDERR "No such file: $fitfile\n";
	    print STDERR "Warning: No fitness file for $org -- do you really want to include it?\n";
	}
	my $seqfile2 = $seqfile . "2";
	AddOrgName($seqfile,$org,$seqfile2) unless NewerThan($seqfile2, $seqfile);
	unless (NewerThan("$seqfile2.pin",$seqfile2)) {
	    system("$blastdir/formatdb -o T -p T -i $seqfile2") == 0 || die "formatdb failed";
	}
	push @seq2files, $seqfile2;
    }

    # Run all pairs of BLASTp jobs
    my @cmds;
    my $jobdir = "$gdir/blast_results";
    mkdir($jobdir) unless  -d $jobdir;
    foreach my $org1 (@orgs) {
	my $db1 = "$gdir/$org1/aaseq2";
	foreach my $org2 (@orgs) {
	    my $db2 = "$gdir/$org2/aaseq2";
	    my $pairFile = "$jobdir/blastp_${org1}_vs_${org2}";
	    unless (NewerThan($pairFile, $db1) && NewerThan($pairFile, $db2)) {
		unlink($pairFile);
		push @cmds, qq{$blastdir/blastall -p blastp -e 1e-5 -z 1e8 -F "m S" -m 8 -i $db1 -d $db2 -o $pairFile.tmp && mv $pairFile.tmp $pairFile};
	    }
	}
    }

    print STDERR "Running BLAST jobs for " . scalar(@cmds) . " pairs of genomes\n";
    if (@cmds > 0) {
	print STDERR "Using $nCPU CPUs\n" if defined $nCPU;
	system("rm -f $jobdir/blastcmds*");
	my $cmdfile = "$jobdir/blastcmds";
	open(CMD, ">", $cmdfile) || die "Cannot write to $cmdfile";
	foreach my $cmd (@cmds) {
	    print CMD "$cmd\n";
	}
	close(CMD) || die "Error writing to $cmdfile";
	my $CPUspec = defined $nCPU ? "-n $nCPU" : "";
	my $submitter = "$Bin/submitter.pl $CPUspec $cmdfile";
	system($submitter);
    }

    # Verify that all files exist and combine them
    my $combfile = "$outdir/aaseqs.blastp";
    open(OUT, ">", $combfile) || die "Cannot write to $combfile";
    foreach my $org1 (@orgs) {
	foreach my $org2 (@orgs) {
	    my $pairFile = "$jobdir/blastp_${org1}_vs_${org2}";
	    open(IN, "<", $pairFile) || die "Cannot read $pairFile -- did BLASTp fail?";
	    while(my $line = <IN>) {
		print OUT $line;
	    }
	    close(IN) || die "Error reading $pairFile";
	}
    }
    close(OUT) || die "Error writing to $combfile";

    print STDERR "Wrote $combfile\n";
    system("$Bin/bbh.pl -blastp $combfile -out $outdir/aaseqs") == 0 || die "bbh.pl failed";
    system("$Bin/paralogs.pl $combfile > $outdir/aaseqs.para") == 0 || die "paralogs.pl failed";
    unlink($combfile);
    print STDERR "Wrote aaseqs.bbh aaseqs.scores aaseqs.para\n";
}

sub AddOrgName($$$) {
    my ($file,$prefix,$file2) = @_;
    open(IN, "<", $file) || die "Error reading $file";
    open(OUT, ">", $file2) || die "Cannot write to $file2";
    while(my $line = <IN>) {
	chomp $line;
        if ($line =~ m/>(.*)$/) {
	    $line = ">$prefix:$1";
	}
	print OUT $line."\n";
    }
    close(IN) || die "Error reading $file";
    close(OUT) || die "Error writing to $file2";
    print STDERR "Wrote $file2\n";
}
