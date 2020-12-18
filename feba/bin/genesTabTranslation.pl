#!/usr/bin/perl -w
use strict;
use Bio::Perl;

die "Usage: genesTabTranslation.pl genes.tab genome.fna > aaseq"
    unless @ARGV == 2;

{
    my ($genesfile, $fnafile) = @ARGV;

    my %seqs = ();
    my $lastName = undef;
    open (FNA, "<", $fnafile) || die "Error reading $fnafile";
    while(<FNA>) {
	s/[\r\n]+$//;
	if (m/>(\S+)/) {
	    $lastName = $1;
	    $seqs{$lastName} = "";
	} else {
	    $seqs{$lastName} .= $_;
	}
    }
    close(FNA) || die "Error reading $fnafile";

    open(GENES, "<", $genesfile) || die "Error reading $genesfile";
    my @header = (); # list of column names
    my %header = (); # column name to index
    my $iScaffold;
    while(<GENES>) {
	chomp;
	my @F = split /\t/, $_, -1; # leave blank columns at end
	if (scalar(@header) == 0) {
	    # read header line
	    @header = @F;
	    %header = map { $header[$_] => $_ } (0..(scalar(@header)-1));	    
	    $iScaffold = $header{"scaffold"};
	    $iScaffold = $header{"scaffoldId"} if !defined $iScaffold;
	    die "No scaffold or scaffoldId column: $_" unless defined $iScaffold;
	    foreach my $col (qw{begin end strand locusId}) {
		die "No $col column" unless exists $header{$col};
	    }
	    if (!exists $header{"type"}) {
		print STDERR "Warning: no type column, assuming all protein-coding\n";
	    }
	} else {
	    # data row
	    my $sc = $F[$iScaffold];
	    my $begin = $F[ $header{"begin"} ];
	    my $end = $F[ $header{"end"} ];
	    my $strand = $F[ $header{"strand"} ];
	    my $locusId = $F[ $header{"locusId"} ];
	    my $type = 1;
	    $type = $F[ $header{"type"} ] if exists $header{"type"};
	    die "Missing fields in $_"
		unless defined $sc && defined $begin && defined $end && defined $strand && defined $locusId && defined $type;
	    die "Unknown scaffold $sc" unless exists $seqs{$sc};
	    die "Invalid begin" unless $begin =~ m/^\d+$/;
	    die "Invalid end" unless $end =~ m/^\d+$/;
	    die "Invalid strand" unless $strand eq "+" || $strand eq "-";
	    next unless $type eq "1";
	    my $seq = $seqs{$sc};
	    my $sclen = length($seq);
	    die "Sequence for $sc is of length $sclen but gene $locusId ends at $end" if $end > $sclen;
	    my $subseq;
	    if ($begin < $end) {
		$subseq = substr($seq, $begin-1, $end-$begin+1);
	    } else {
		print STDERR "Warning: wraparound orf $sc $begin $end\n";
		$subseq = substr($seq, $begin-1) . substr($seq, 0, $end);
	    }
	    if ($strand eq "-") {
		$subseq = reverse_complement_as_string($subseq);
	    }
	    my $len = length($subseq);
	    if (($len % 3) != 0) {
		print STDERR "Skipping gene of irregular length $locusId at $sc $begin:$end $strand";
	    }
	    my $aaseq = translate_as_string($subseq);
	    $aaseq =~ s/[*]$//; # ignore stop codon at end, if present
	    print STDERR "Warning: stop codon within $locusId\n" if $aaseq =~ m/[*]/;
	    print ">$locusId\n$aaseq\n";
	}
    }
    close(GENES) || die "Error reading $genesfile";
}
