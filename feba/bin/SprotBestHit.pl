#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for NewerThan

my $sprot = "sprot.faa";
my $gdir = "g";
my $hitdir = "g/blast_results";
my $minIdentity = 30;
my $minCoverage = 0.80;
my $outdir = ".";
my $usage = <<END
Usage: SprotBestHit.pl organism1 ... organismN > besthit.sprot
Optional arguments:
   -out directory to write the output files in
   -sprot sprot.faa -- the SwissProt (curated UniProt) database in fasta format
   -gdir $gdir
   -hitdir $hitdir -- where to look for the rapsearch2 results
   -minIdentity $minIdentity -- minimum %identity
   -minCoverage $minCoverage -- minimum coverage both ways

Assumes that hitdir/sprot_org.m8 contains the rapsearch2 results for that
aaseq2 file versus swissprot, from a command like

rapsearch -q g/org/aaseq2 -t a -e -2 -b 30 -v 30 -d sprot.rap -o g/blast_results/sprot_org

Writes two tab-delimited files:
db.BestHitSwissProt has fields
orgId, locusId, swissprot_accession, swissprot_id, identity

db.SwissProtDesc has fields
swissprot_accession, swissprot_id, geneName (or empty), desc, organism
END
    ;

die $usage unless GetOptions('sprot=s' => \$sprot,
                             'gdir=s' => \$gdir,
                             'hitdir=s' => \$hitdir,
                             'minIdentity=f' => \$minIdentity,
                             'minCoverage=f' => \$minCoverage,
                             'out=s' => \$outdir );
die "Must specify which organisms to process:\n$usage" if @ARGV == 0;;
die "No such directory: $outdir\n" unless -d $outdir;
my @orgs = @ARGV;
my $maxLog10E = -2;

# Check for inputs
die "No such file: $sprot" unless -e $sprot;

foreach my $org (@orgs) {
    my $hitfile = "$hitdir/sprot_$org.m8";
    die "No SwissProt results for $org, no $hitfile" unless -e $hitfile;
    my $aafile = "$gdir/$org/aaseq2";
    die "No aaseq2 for $org\n" unless -e $aafile;
    die "SwissProt results for $org are out of date\n" unless NewerThan($hitfile,$aafile);
}

# Parse the swissprot a.a. file to save the description part of the header
# and the length of each sequence.

# The key is the fasta id, i.e. "sp|Q6GZW9|006R_FRG3G"
my %sprotLen = ();

open(SPROT, "<", $sprot) || die "Cannot read $sprot";
my $name = undef;
while(<SPROT>) {
    chomp;
    if (m/>(.*)/) {
        my $header = $1;
        $header =~ m/^([0-9A-Za-z_|]+) +(.*)$/ || die "Cannot parse SwissProt description line:\n$header\n ";
        $name = $1;
        $sprotLen{$name} = 0;
    } else {
        $sprotLen{$name} += length($_);
    }
}
close(SPROT) || die "Error reading $sprot";
print STDERR "Parsed " . scalar(keys %sprotLen) . " sequences from $sprot\n";

my %sprotUsed = ();

# Identify best hits and print to db.BestHitSwissProt
my $bhFile = "$outdir/db.BestHitSwissProt";
open(BESTHIT, ">", $bhFile) || die "Cannot write to $bhFile";
foreach my $org (@orgs) {
    print STDERR "Processing $org\n";
    my $hitfile = "$hitdir/sprot_$org.m8";
    my $aafile = "$gdir/$org/aaseq2";

    # First, the length of each sequence
    my %qlen = ();
    open(FAA, "<", $aafile) || die "Cannot read $aafile";
    my $seqname = undef;
    while(<FAA>) {
        chomp;
        if (m/>(.*)/) {
            $seqname = $1;
            $qlen{$seqname} = 0;
        } else {
            $qlen{$seqname} += length($_);
        }
    }
    close(FAA) || die "Error reading $aafile";
    while (my ($seq,$len) = each %qlen) {
        die "No sequence for $seq" unless $len > 0;
    }
    print STDERR "Read " . scalar(keys %qlen) . " sequences for $org\n";

    my %besthit = (); # query => [id, bits, identity]
    open(HITS, "<", $hitfile) || die "Cannot read $hitfile";
    while(my $line = <HITS>) {
        next if $line =~ m/^#/; # header line
        chomp $line;
        my ($query,$subject,$identity,$alnlen,$mismatch,$gaps,$qstart,$qend,$sstart,$send,$log10E,$bits) = split /\t/, $line;
        die "Cannot parse\n$_\nfrom $hitfile" unless defined $bits;
        next unless $identity >= $minIdentity && $log10E <= $maxLog10E;
        next if exists $besthit{$query} && $besthit{$query}->[1] >= $bits;
        die "No length for $query" unless exists $qlen{$query};
        die "Invalid length for $query" if $qend > $qlen{$query};
        die "No length for $subject" unless exists $sprotLen{$subject};
        die "Invalid length for $subject, $send vs. $sprotLen{$subject}" if $send > $sprotLen{$subject};
        if (($qend-$qstart+1) >= $minCoverage * $qlen{$query}
            && ($send-$sstart+1) >= $minCoverage * $sprotLen{$subject}) {
            $besthit{$query} = [$subject,$bits,$identity];
        }
    }
    close(HITS) || die "Error reading $hitfile";
    print STDERR "Found good hits for " . scalar(keys %besthit) . " of " . scalar(keys %qlen) . " proteins in $org\n";
    while (my ($query, $row) = each %besthit) {
        my ($subject,undef,$identity) = @$row;
        my ($orgId,$locusId) = split /:/, $query;
        die "Cannot parse query $query" unless defined $locusId;
        my ($base,$sp_acc,$sp_id) = split /[|]/, $subject;
        die "Cannot parse subject $subject" unless defined $sp_id;
        print BESTHIT join("\t", $orgId, $locusId, $sp_acc, $sp_id, $identity)."\n";
        $sprotUsed{$subject} = 1;
    }
}
close(BESTHIT) || die "Error writing to $bhFile";
print STDERR "Processed " . scalar(@orgs) . " organisms to $bhFile\n";

open(SPROT, "<", $sprot) || die "Cannot read $sprot";
my $descfile = "$outdir/db.SwissProtDesc";
open(DESC, ">", $descfile) || die "Cannot write to $descfile";
while(<SPROT>) {
    next unless m/^>/;
    chomp;
    s/^>//;
    m/^([0-9A-Za-z_|]+) +(.*)$/ || die;
    my $id = $1;
    next unless exists $sprotUsed{$id};
    my $desc = $2;
    $desc =~ m/^([^=]+) OS=(.*)$/ || die "Cannot parse organism name from $desc";
    $desc = $1;
    my $tail = $2;
    my $org = "";
    my $gene = "";
    $tail =~ s/ OX=\d+\b//; # remove the organism identifier
    if ($tail =~ m/^([^=]+) GN=(\S+)/ ) {
        $org = $1;
        $gene = $2;
    } elsif ($tail =~ m/^([^=]+) PE=/) {
        $org = $1;
    } else {
        die "Cannot parse tail part of sprot desc line $tail";
    }
    my (undef,$sp_acc,$sp_id) = split /[|]/, $id;
    die $id unless defined $sp_id;
    print DESC join("\t", $sp_acc, $sp_id, $gene, $desc, $org)."\n";
}
close(DESC) || die "Error writing to $descfile";
print STDERR "Wrote $descfile\n";
