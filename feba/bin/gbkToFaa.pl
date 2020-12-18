#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

{
    die "Usage: gbkToFaa.pl GeneTable GenbankFile > aaFastaFile\n"
	unless @ARGV==2;
    my ($genes,$file) = @ARGV;

    my $obj = Bio::SeqIO->new(-file => $file, -format => "genbank");

    my %tagToId = (); # sysName + begin => locusId
    open(GENES, "<", $genes) || die "Cannot read $genes";
    my $header = <GENES>;
    chomp $header;
    my @header = split /\t/, $header;
    my %headerNames = map { $header[$_] => $_ } (0..(scalar(@header)-1));
    die "Missing information for sysName" unless exists $headerNames{'sysName'};
    die "Missing information for locusId" unless exists $headerNames{'locusId'};
    die "Missing information for locusId" unless exists $headerNames{'begin'};
    die "Missing information for locusId" unless exists $headerNames{'type'};

    while(my $line = <GENES>) {
	chomp $line;
	my @F = split /\t/, $line, -1;
	my $locusId = $F[ $headerNames{'locusId'} ];
	my $sysName = $F[ $headerNames{'sysName'} ];
	my $begin = $F[ $headerNames{'begin'} ];
	my $type = $F[ $headerNames{'type'} ];
	next unless defined $sysName && defined $locusId && $sysName ne "";
        next unless $type eq "1";
	my $tag = $sysName . ":" . $begin;
	die "Duplicate entry for $tag" if exists $tagToId{$tag};
	$tagToId{$tag} = $locusId;
    }
    close(GENES) || die "Error reading from $genes";

    my %locusTags = (); # for tracking duplicates
    while(my $seq = $obj->next_seq) {
	foreach my $feature ($seq->get_SeqFeatures) {
	    if ($feature->primary_tag eq "CDS") {
		my $sysName = ($feature->get_tag_values("locus_tag"))[0];
		my $start = $feature->start;
		my $end = $feature->end;
		my $begin = $start < $end ? $start : $end;
		my $tag = $sysName . ":" . $begin;
		print "Warning: duplicate locus $tag" if exists $locusTags{$sysName};
		$locusTags{$tag} = 1;
		my $locusId = $tagToId{$tag};
		if (!defined $locusId) {
		    print STDERR "No gene entry for locustag:begin $tag (OK if wrap-around ORF or pseudo)\n";
		} else {
                    my $aaseq = ($feature->get_tag_values("translation"))[0];
		    print ">$locusId\n$aaseq\n";
		}
	    }
	}
    }
    print STDERR "Wrote a.a. sequences for " . scalar(keys %locusTags) . " genes\n";
}
