#!/usr/bin/perl -w
use strict;
# Convert a gff file from JGI or RefSeq to a table of genes
# In RefSeq files, the locus_tag is on a gene line that is just before the CDS/tRNA/rRNA line

{
    my $prefix = "GFF";
    if (@ARGV >= 2 && $ARGV[0] eq "-prefix") {
	shift @ARGV;
	$prefix = shift @ARGV;
    }
    die "gffToGenes.pl [ -prefix GFF ] < file.gff > out.gff\n" unless @ARGV == 0;
    my $nGFF = 0;
    my %types = ("CDS" => 1, "rRNA" => 2, "tRNA" => 5, "RNA" => 6, "transcript" => 6);
    print join("\t",qw{locusId sysName type scaffoldId begin end strand name desc})."\n";
    my %sysNameSeen = ();
    my $last = undef; # the preceding gene line, to get locus_tag from
    while(<STDIN>) {
	chomp;
	next if m/^#[#!]/;
	my ($scaffold, $source, $type, $begin, $end, $score, $phase, undef, $attributes) = split /\t/, $_;
	if ($type eq "gene") {
	    $last = { 'begin' => $begin, 'end' => $end, 'attributes' => $attributes };
	}
	next unless exists $types{$type}; # skip gene, exon, etc. lines
	die "Cannot parse phase $phase" unless $phase =~ m/^-?1$/ || $phase eq "+" || $phase eq "-";
	unless ($begin < $end) {
	    print STDERR "Warning: skipping GFF entry because begin >= end: "
                . join(" ", $scaffold, $begin, $end, $phase, $attributes)."\n";
	    next;
	}
        if (abs($begin - $end) > 50000) {
            print STDERR "Warning: skipping huge GFF entry (possibly a wraparound ORF): "
                . join(" ", $scaffold, $begin, $end, $phase, $attributes)."\n";
            next;
        }
	$nGFF++;
	my $strand = $phase;
	$strand = "+"if $phase eq "1";
	$strand = "-"if $phase eq "-1";
	my $desc = "";
	if ($attributes =~ m/product='(.*)'$/) { # IMG does not quote ' characters within product, use end of line
	    $desc = $1;
	} elsif ($attributes =~ m/product=(\S+)$/) { # IMG does not quote short descriptions at all; hope they are at end of line
	    $desc = $1;
	} elsif ($attributes =~ m/product='([^']+)' /) { # allow product not at end if quoted and contains no quotes
	    $desc = $1;
	} elsif ($attributes =~ m/ product=([^ ]+) [a-zA-Z_]+=/) {
	    $desc = $1;
	} elsif ($attributes =~ m/ product='([^=]+)' protein_id=/) {
	    $desc = $1; # a hack for improperly quoted products with protein_id afterwards
	} elsif ($attributes =~ m/product=(.*)$/) {
	    $desc = $1;
	} else {
	    print STDERR "Cannot parse product for $type from: $attributes\n";
	    $desc = $type;
	}
	$desc =~ s/;[^ ].*$//; # in case parsing above fails to separate out the fields
	# deal with occasional use of HTML entities
	$desc =~ s/%2C/ /g;
	$desc =~ s/%3B/;/g;

	my $sysName = "";
	# \S for non-whitespace, \b for word boundary
	$sysName = $1 if $attributes =~ m/\blocus_tag=([^ ;]+)\b/;
	if ($sysName eq "" && defined $last && $last->{begin} == $begin && $last->{end} == $end) {
	    $sysName = $1 if $last->{attributes} =~ m/\blocus_tag=([^ ;]+)\b/;
	}

	my $locusId = "$prefix$nGFF";
	if ($sysName ne "") {
	    if (exists $sysNameSeen{$sysName}) {
		print STDERR "Warning: duplicate locus tag $sysName\n";
	    } else {
		$locusId = $sysName;
	    }
	    $sysNameSeen{$sysName} = 1;
	}

	# output types: 1 for protein, 2 for rRNA, 5 for tRNA, 6 for other non-coding RNA,
	# 7 for pseudogene of protein-coding genes, 8 for pseudogene of RNA gene
	my $out_type = $types{$type};
	my $isPseudo = $attributes =~ m/pseudo=true/
            || $attributes =~ m/partial=true/
            || $attributes =~ m/pseudo=_no_value/;
	$out_type = 7 if $isPseudo && $out_type == 1;
	$out_type = 8 if $isPseudo && $out_type == 8;
	# a hack to identify (most of) the JGI pseudogenes, as they are not indicated in the gff file.
	if ($out_type == 1 && ($end-$begin+1) % 3 != 0) {
	    print STDERR "Warning: extent $begin:$end not a multiple of 3 for CDS $locusId $sysName => pseudogene\n";
	    $out_type = 7;
	}
	print join("\t", $locusId, $sysName, $out_type, $scaffold, $begin, $end, $strand, "", $desc)."\n";
    }
}
