#!/usr/bin/perl -w
# Given the output of MapTnSeq.pl, possibly from multiple runs, choose the reliable and unambiguous tags
# Makes a tab-delimited table
#
# Input is assumed to contain
# read, barcode, scaffold, position, strand, unique flag, qBeg, qEnd, bit score, %identity
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for reverseComplement
sub Variants($); # return all 1-nt variants for a sequence

my $minN = 10;
my $minFrac = 0.75;
my $minRatio = 8.0;
my $maxQBeg = 3;

my $usage = <<END
Usage: DesignRandomPool.pl -pool pool_file -genes genes_table [ -minN $minN ]
          [ -minFrac $minFrac ] [ -minRatio $minRatio ] [ -maxQBeg $maxQBeg ]
          MapTnSeq_file1 ... MapTnSeq_fileN

DesignRandomPool.pl identifies the reliably mapped barcodes and writes
a pool file, as well as to auxilliary files pool.hit (strains per
gene), pool.unhit (proteins without insertions), pool.surprise.

The MapTnSeq files must be tab-delimited with fields
    read,barcode,scaffold,pos,strand,uniq,qBeg,qEnd,score,identity
where qBeg and qEnd are the positions in the read, after the
transposon sequence is removed, that match the genome.
Gzipped (*.gz) input files are also supported.

The genes table must include the fields
   scaffoldId, begin, end, strand, desc
and it should either include only protein-coding genes or it should
include the field 'type' with type=1 for protein-coding genes.

Optional arguments:

minN is the minimum number of "good" reads for a barcode supporting
its mapping.  (Good means a unique hit to the genome and qBeg=1.)

minFrac is the minimum fraction of input reads for the barcode that
agree with the preferred mapping.

minRatio is the minimum ratio of
reads for preferred mapping over the 2nd-most-frequent mapping.
END
    ;

{
    my $poolfile = undef;
    my $genesfile = undef;
    GetOptions('minN=i' => \$minN,
               'minFrac=f' => \$minFrac, 'minRatio=f' => \$minRatio,
	       'pool=s' => \$poolfile,
	       'genes=s' => \$genesfile,
               'maxQBeg=i' => \$maxQBeg )
        || die $usage;
    die $usage if !defined $poolfile;
    die $usage if !defined $genesfile;
    die "No such file: $genesfile" unless -e $genesfile;
    die "minN must be at least 2\n" unless $minN >= 2;

    my $nMapped = 0; # reads considered
    my $nSkipQBeg = 0; # reads ignored because qBeg > maxQBeg

    my %barPosCount = (); # barcode => scaffold:strand:position => list of nReads, nGoodReads
    # where a "good" read has uniq=1 and qBeg=1

    my %pastEnd = (); # number of reads for that barcode mapped past the end
    # (pastEnd reads are not included in barPosCount)

    print STDERR join("\n", "Reading mapping files:", @ARGV)."\n"
      if @ARGV > 0;

    foreach my $file (@ARGV) {
      my $fhIn;
      if ($file =~ m/[.]gz$/) {
        open($fhIn, '-|', 'zcat', $file) || die "Cannot run zcat on $file";
      } else {
        open($fhIn, "<", $file) || die "Cannot read $file";
      }
      while(<$fhIn>) {
        chomp;
        my ($read,$barcode,$scaffold,$pos,$strand,$uniq,$qBeg,$qEnd,$score,$identity) = split /\t/, $_;
        if ($scaffold eq "pastEnd") {
            $pastEnd{$barcode}++;
            $nMapped++;
        } elsif ($maxQBeg >= 1 && $qBeg <= $maxQBeg)  {
	    my $key = join(":",$scaffold,$strand,$pos);
	    $barPosCount{$barcode}{$key}[0]++;
	    $barPosCount{$barcode}{$key}[1]++ if $uniq eq "1" && $qBeg eq "1";
            $nMapped++;
        } else {
            $nSkipQBeg++;
        }
      }
      close($fhIn) || die "Error reading $file";
    }
    print STDERR "Read $nMapped mapped reads for " . scalar(keys %barPosCount) . " distinct barcodes\n";
    print STDERR "(Skipped $nSkipQBeg reads with qBeg > $maxQBeg)\n" if $nSkipQBeg > 0;

    open(POOL, ">", $poolfile) || die "Cannot write to $poolfile";
    print POOL join("\t", "barcode", "rcbarcode", "nTot",
               "n","scaffold","strand","pos",
               "n2","scaffold2","strand2","pos2","nPastEnd")."\n";

    # The Chao2 estimator of diversity is tot + f1**2/(2*f2), where tot is the number of barcodes seen,
    # f1 is the number seen exactly once, and f2 is the number seen exactly twice
    # This does not account for sequencing error and does not account for reuse of a barcode across strains
    my $f1 = 0;
    my $f2 = 0;
    my $totcodes = 0;

    my %nInCategory = ("Usable" => 0); # classification of barcodes
    my ($SCAFFOLD,$POS,$STRAND,$UNIQ,$QBEG,$QEND,$SCORE,$IDENTITY) = (0,1,2,3,4,5,6,7);
    my $nMulti = 0; # barcodes seen >once
    my %barcodeAt = (); # barcode to list of nTot, nMax, at, nNext, nextAt
    my $nReadsForUsable = 0;
    while(my ($barcode, $hash) = each %barPosCount) {
        my $nPastEnd = $pastEnd{$barcode} || 0;
	my $nTot = $nPastEnd;
	while (my ($key,$value) = each %$hash) {
	    $nTot += $value->[0]; # the total reads
	}
	if ($nTot == 1) {
	    $f1++;
	} elsif ($nTot == 2) {
	    $f2++;
	}
	$totcodes++;
        next unless $nTot >= $minN;
        $nMulti++;

        # Is there a location that accounts for most of the reads, including the past-end reads?
        my @at = sort { $hash->{$b}[0] <=> $hash->{$a}[0] } keys(%$hash);
        my $nMax = @at == 0 ? 0 : $hash->{$at[0]}[0];
        if ($nPastEnd >= $nTot/2 || $nPastEnd >= $nMax) {
            $nInCategory{"PastEnd"}++;
            my $n2 = $nMax || 0; # note we do not report secondary location (doubt it matters)
            print POOL join("\t", $barcode, &reverseComplement($barcode),
                       $nTot, $nPastEnd, "pastEnd","","",$n2,"","","",$nPastEnd)."\n";
            next;
        }
        unless($nMax >= $minN && $nMax/$nTot >= $minFrac) {
            $nInCategory{"NoCons"}++;
            next;
        }

        my $maxAt = $at[0];

        # checking unique & qbeg=1 -- but the latter part may be redundant with maxQBeg
        my $nGood = $hash->{$maxAt}[1] || 0;
        unless ($nGood >= $minN && $nGood/$nTot >= $minFrac) {
            $nInCategory{"FewGood"}++;
            next;
        }

        my $nextAt = "::";
        my $nNext = 0;
        if (@at > 1) {
          $nextAt = $at[1];
          $nNext = $hash->{$nextAt}[0];
        }
        unless($nMax >= $minRatio * $nNext) {
            $nInCategory{"LoRatio"}++;
            next;
        }
        $barcodeAt{$barcode} = [ $nTot, $nMax, $maxAt, $nNext, $nextAt ];
        $nInCategory{"Usable"}++;
        $nReadsForUsable += $nTot;
    }

    print STDERR "$nMulti barcodes seen $minN or more times, map $nInCategory{Usable} (minFrac $minFrac minRatio $minRatio)\n";
    foreach my $category (sort keys %nInCategory) {
        print STDERR sprintf("%s\t%d\t%.4f\n", $category, $nInCategory{$category},
			     $nInCategory{$category} / $nMulti)
	    if $nInCategory{$category} > 0;
    }

    my $nOut = 0;
    my $nMasked = 0;
    my $nMaskedReads = 0;
    while (my ($barcode,$row) = each %barcodeAt) {
        my ($nTot,$nMax,$maxAt,$nNext,$nextAt) = @$row;
        my @variants = Variants($barcode);
        my $mask = 0;
        foreach my $variant (@variants) {
            if (exists $barcodeAt{$variant}
                && $barcodeAt{$variant}[0] > $nTot
                && $barcodeAt{$variant}[2] eq $maxAt) {
                $nMasked++;
                $nMaskedReads += $nMax;
                $mask = 1;
                next;
                last;
            }
        }
        next if $mask;
        my @atSplit = split /:/, $maxAt;
        my @nextAtSplit = split /:/, $nextAt, -1; # keep empty entries
        my $nPastEnd = $pastEnd{$barcode} || 0;
        print POOL join("\t", $barcode, reverseComplement($barcode),
                   $nTot, $nMax, @atSplit, $nNext, @nextAtSplit, $nPastEnd)."\n";
        $nOut++;
    }
    close(POOL) || die "Error writing to $poolfile";
    print STDERR "Masked $nMasked off-by-1 barcodes ($nMaskedReads reads) leaving $nOut barcodes\n";
    print STDERR sprintf("Reads for those barcodes: %d of %d (%.1f%%)\n",
                         $nReadsForUsable, $nMapped, 100*$nReadsForUsable/($nMapped + 1e-6));
    my $chao = int($totcodes + $f1**2/(2*$f2 + 1));
    print STDERR "Chao2 estimate of #barcodes present (may be inflated for sequencing error): $chao\n";
    system("$Bin/../lib/PoolStats.R", $poolfile, $genesfile, $nMapped) == 0 || die $!;
}

sub Variants($) {
    my ($baseseq) = @_;
    my @out = ();
    $baseseq = uc($baseseq);
    my $len = length($baseseq);
    foreach my $i (0..($len-1)) {
        my $pre = substr($baseseq,0,$i);
        my $char = substr($baseseq,$i,1);
        my $post = substr($baseseq,$i+1);
        next unless $char eq "A" || $char eq "C" || $char eq "G" || $char eq "T";
        foreach my $newchar (qw{A C G T}) {
            push @out, $pre . $newchar . $post unless $newchar eq $char;
        }
    }
    return(@out);
}

# no longer used
sub uniq(@) {
    my %seen = ();
    foreach my $value (@_) { $seen{$value} = 1; }
    return keys(%seen);
}
