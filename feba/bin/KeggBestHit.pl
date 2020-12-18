#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for NewerThan

my $kegg = "kegg.faa";
my $gdir = "g";
my $hitdir = "g/blast_results";
my $minIdentity = 30;
my $minCoverage = 0.80;
my $usage = <<END
Usage: KeggBestHit.pl organism1 ... organismN > besthit.kegg
Optional arguments:
   -kegg kegg.faa -- the KEGG database in fasta format
       (This script works with the last public release from June 2011)
   -gdir $gdir
   -hitdir $hitdir -- where to look for the rapsearch2 results
   -minIdentity $minIdentity -- minimum %identity
   -minCoverage $minCoverage -- minimum coverage both ways

Assumes that hitdir/KEGG_org.m8 contains the rapsearch2 results for that
aaseq2 file versus KEGG, from a command like

rapsearch -q g/org/aaseq2 -t a -e -2 -b 30 -v 30 -d kegg.rap -o g/blast_results/KEGG_org

Writes a tab-delimited file (no header) with fields
orgId:locusId, keggId, identity, K group numbers, K group descriptions

where K group numbers (if more than 1) are joined by ","
and K group descriptions are joined by "; "
END
    ;

die $usage unless GetOptions('kegg=s' => \$kegg,
                             'gdir=s' => \$gdir,
                             'hitdir=s' => \$hitdir,
                             'minIdentity=f' => \$minIdentity,
                             'minCoverage=f' => \$minCoverage);
die "Must specify which organisms to process:\n$usage" if @ARGV == 0;;
my @orgs = @ARGV;
my $maxLog10E = -2;

my %kegglen = (); # id => #a.a.
my %keggdesc = (); # id => description
my %keggko = (); # id => list of [KO id, KO description]

my %favoredOrg = ("eco" => 1); # prefer hits to these KEGG organism ids

# Check for inputs
die "No such file: $kegg" unless -e $kegg;

foreach my $org (@orgs) {
    my $hitfile = "$hitdir/KEGG_$org.m8";
    die "No KEGG results for $org, no $hitfile" unless -e $hitfile;
    my $aafile = "$gdir/$org/aaseq2";
    die "No aaseq2 for $org\n" unless -e $aafile;
    die "KEGG results for $org are out of date\n" unless NewerThan($hitfile,$aafile);
}

# Parse KEGG
# Note, "No updates are made to original data, such as gene names and
# descriptions given by RefSeq or GenBank, even if they are
# inconsistent with the KO assignment."
# So I will use the K descriptions only.

open(KEGG, "<", $kegg) || die "Cannot read $kegg";
my $name = undef;
while(<KEGG>) {
    chomp;
    if (m/>(.*)/) {
        my $header = $1;
        $header =~ m/(\S+:\S+) +(.*)$/ || die "Cannot parse KEGG description line:\n$header\n ";
        $name = $1;
        my $descs = $2;
        $kegglen{$name} = 0;
        $keggko{$name} = [];
        my @desc = split /; /, $descs;
        die "Cannot parse descs for $name:\n$descs\n " unless @desc >= 1;
        # Usually the first one is the description; the others are ortholog groups, there can be more than one of those
        # Or, the first one is a gene name, the next is a description
        $keggdesc{$name} = shift @desc;
        foreach my $desc (@desc) {
            if ($desc =~ m/^(K\d+) (.*)$/) { #K number and description
                push @{ $keggko{$name} }, [ $1, $2 ];
            } elsif ($desc =~ m/^(K\d+)$/) { #K number only
                push @{ $keggko{$name} }, [ $1, "" ];
            } else {
                # not a K number, is an additional annotation instead
                $keggdesc{$name} = $keggdesc{$name} . "; " . $desc;
            }
        }
    } else {
        $kegglen{$name} += length($_);
    }
}
close(KEGG) || die "Error reading $kegg";
print STDERR "Parsed " . scalar(keys %kegglen) . " sequences from $kegg\n";

foreach my $org (@orgs) {
    print STDERR "Processing $org\n";
    my $hitfile = "$hitdir/KEGG_$org.m8";
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

    my %besthit = (); # query => [KEGG id, bits, identity]
    open(HITS, "<", $hitfile) || die "Cannot read $hitfile";
    while(my $line = <HITS>) {
        next if $line =~ m/^#/; # header line
        chomp $line;
        my ($query,$subject,$identity,$alnlen,$mismatch,$gaps,$qstart,$qend,$sstart,$send,$log10E,$bits) = split /\t/, $line;
        die "Cannot parse\n$_\nfrom $hitfile" unless defined $bits;
        next unless $identity >= $minIdentity && $log10E <= $maxLog10E;
        my $old = $besthit{$query};
        if (defined $old) {
            if ($identity > 90) {
                # try to prefer hits to organisms in favored org if they are not noticeably worse
                next if $old->[1] > $bits + 4; # old hit is noticeably better
                my $oldsub = $old->[0];
                my $oldorg = $oldsub; $oldorg =~ s/:.*//;
                next if exists $favoredOrg{$oldorg}; # favoring this hit
            } else {
                next if $old->[1] >= $bits; # already have a better hit
            }
        }
        die "No length for $query" unless exists $qlen{$query};
        die "Invalid length for $query" if $qend > $qlen{$query};
        die "No length for $subject" unless exists $kegglen{$subject};
        die "Invalid length for $subject, $send vs. $kegglen{$subject}" if $send > $kegglen{$subject};
        if (($qend-$qstart+1) >= $minCoverage * $qlen{$query}
            && ($send-$sstart+1) >= $minCoverage * $kegglen{$subject}) {
            $besthit{$query} = [$subject,$bits,$identity];
        }
    }
    close(HITS) || die "Error reading $hitfile";
    print STDERR "Found good hits for " . scalar(keys %besthit) . " of " . scalar(keys %qlen) . " proteins in $org\n";
    foreach my $query (sort keys %besthit) {
        my ($subject,$bits,$identity) = @{ $besthit{$query} };
        my @ko = @{ $keggko{$subject} };
        my @knums = map {$_->[0]} @ko;
        my @kdesc = map {$_->[1]} @ko;
        print join("\t", $query, $subject, $identity, join(",", @knums), join("; ", @kdesc))."\n";
    }
}
print STDERR "Processed " . scalar(@orgs) . " organisms\n";
