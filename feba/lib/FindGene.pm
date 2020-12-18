
package FindGene;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( LocationToGene CheckGeneLocations );

# Given scaffold, pos, and hash of scaffold to list of sorted genes,
# returns the locusId and the fraction through the gene it is in (as a list of 2 elements)
#
# If the location overlaps multiple genes or no genes, returns locusId = "", f = "".
#
# Each gene should be a hash that contains begin, end, strand, and locusId
#
# This code does not support orfs that wrap around the origin, and it may not give correct
# results if there are complicated overlaps between ORFs. In particular, it only checks
# the adjacent ORF on either side to see if there is an overlap.
sub LocationToGene($$$) {
    my ($scaffold, $pos, $sortedGenes) = @_;
    return ("","") if $scaffold eq "pastEnd";

    my $genelist = $sortedGenes->{$scaffold};
    return ("","") if !defined $genelist;

    # binary search
    # at all times, either the true index is between lo and hi, or there is no hit
    my $nGenes = scalar(@$genelist);
    my $lo = 0;
    my $hi = $nGenes-1;
    for (my $nRound = 0; $nRound < 100000; $nRound++) {
	my $mid = int(($lo+$hi)/2);
	my $iBegin = $genelist->[$mid]{begin};
	my $iEnd = $genelist->[$mid]{end};

	if ($pos < $iBegin) {
	    return ("","") if $mid == $lo;
	    $hi = $mid-1;
	} elsif ($pos > $iEnd) {
	    return ("","") if $mid == $hi;
	    $lo = $mid+1;
	} else {
	    # does the previous or next gene also overlap this position?
	    return ("","") if ($mid > 0 && $genelist->[$mid-1]{begin} <= $pos && $pos <= $genelist->[$mid-1]{end});
	    return ("","") if ($mid < $nGenes-1 && $genelist->[$mid+1]{begin} <= $pos && $pos <= $genelist->[$mid+1]{end});
	    my $f = $iBegin == $iEnd ? 0 : ($pos - $iBegin)/($iEnd-$iBegin);
	    my $strand = $genelist->[$mid]{strand};
	    # insertions near N terminus of gene should have f near 0 regardless of strand
	    $f = 1.0-$f if $strand eq "-";
	    return($genelist->[$mid]{locusId}, $f);
	}
    }
    die "Unreachable";
}

# Given a hash of scaffold to a list of sorted genes, check to see that
# for almost all genes, the center of the gene does not overlap other genes
# and is identified as the hit for that gene
#
# Each gene should be a hash that contains begin, end, strand, and locusId
#
# Writes a short report to STDERR and reports a reference to the list of locusIds that fail
# (with overlap-adjacent cases excluded)
sub CheckGeneLocations($) {
    my ($sortedGenes) = @_;
    my @ok = ();
    my @wrap = ();
    my @overlap = ();
    my @fail = ();
    while (my ($scaffoldId,$geneList) = each %$sortedGenes) {
	foreach my $i (0..(scalar(@$geneList)-1)) {
	    my $gene = $geneList->[$i];
	    if ($gene->{end} < $gene->{begin}) {
		push @wrap, $gene->{locusId};
	    } else {
		my $pos = int( ($gene->{begin} + $gene->{end} + 1)/2 );
		my ($locusHit, $f) = LocationToGene($scaffoldId, $pos, $sortedGenes);
		if ($locusHit eq "") {
		    # does this position overlap adjacent genes?
		    if ($i > 0 &&  $geneList->[$i-1]{end} >= $pos) {
			push @overlap, $gene->{locusId};
		    } elsif ($i < scalar(@$geneList)-1 && $geneList->[$i+1]{begin} <= $pos) {
			push @overlap, $gene->{locusId};
		    } else {
			push @fail, $gene->{locusId};
		    }
		} elsif ($locusHit eq $gene->{locusId}) {
		    push @ok, $gene->{locusId};
		} else {
		    # print STDERR sprintf("Gene %s hits %s instead!\n", $gene->{locusId}, $locusHit);
		    push @fail, $gene->{locusId};
		}
	    }
	}
    }
    print STDERR sprintf("Check gene locations: success %d wrap %d overlap-adjacent %d failure %d\n",
			 scalar(@ok), scalar(@wrap), scalar(@overlap), scalar(@fail));
    print STDERR "Failures: " . join(" ", @fail)."\n"
	if scalar(@fail) > 0;
    return \@fail;
}

1;
