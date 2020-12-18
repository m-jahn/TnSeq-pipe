#!/usr/bin/perl -w
use strict;

# paralogs.pl converts the BLASTp results, with sequence identifiers of the form
# organism:locusId, into a list of paralogous pairs, of the form
# locus1, locus2, bits, ratio, where ratio = bits/selfscore(locus1).
# Hits are reported in whichever direction(s) they appear in the BLASTp results.
# Only the best hit for each pair is reported

my $minRatio = 0.25; # ignore hits beneath this ratio

{
    die "Usage: paralogs.pl aaseqs.blastp > aaseqs.para\n" unless @ARGV==1;
    my $infile = $ARGV[0];

    # First pass: compute self scores
    my %selfScore = ();
    open(IN, "<", $infile) || die "Cannot read $infile";
    while(<IN>) {
	chomp;
	my ($id1,$id2,$identity,$match,$mm,$gap,$qbeg,$qend,$sbeg,$send,$eval,$bits) = split /\t/, $_;
	die "Wrong number of columns in input" unless defined $bits;
	die "Invalid bit score $bits" unless $bits =~ m/^ *[0-9e.+]+$/;
	if ($id1 eq $id2) {
	    $selfScore{$id1} = $bits if !exists $selfScore{$id1} || $selfScore{$id1} < $bits;
	}
    }
    close(IN) || die "Error reading $infile";

    # Second pass: report paralogs
    open(IN, "<", $infile) || die "Cannot read $infile";
    while(<IN>) {
	chomp;
	my ($id1,$id2,$identity,$match,$mm,$gap,$qbeg,$qend,$sbeg,$send,$eval,$bits) = split /\t/, $_;
	my ($org1,$locus1) = split /:/, $id1;
	my ($org2,$locus2) = split /:/, $id2;
	if ($org1 eq $org2 && $locus1 ne $locus2) {
	    if (!exists $selfScore{$id1}) {
		print STDERR "Warning: $id1 to $id2 score $bits but no hit $id1 to itself\n";
	    } else {
		my $ratio = $bits / $selfScore{$id1};
		print join("\t", $org1, $locus1, $locus2, $bits, $ratio)."\n" if $ratio >= $minRatio;
	    }
	}
    }
    close(IN) || die "Error reading $infile";
}
