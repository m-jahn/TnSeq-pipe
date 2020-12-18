#!/usr/bin/perl -w
use strict;

# This script handles genbank files that have both CONTIG and ORIGIN sections
# that fail with bioperl's parser (which is used by gbkToSeq.pl)

{
    die "Usage: gbkToSeq2.pl GenbankFile > fastaFile\n"
	unless @ARGV==1;
    my ($file) = @ARGV;
    open(FILE, "<", $file) || die "Cannot read $file\n";
    my $locus = undef;
    my $seqsofar = undef;
    while(my $line = <FILE>) {
        chomp $line;
        $line =~ s/ +$//;
        if ($line =~ m/^LOCUS +([A-Za-z0-9._]+) /) {
            die "Second LOCUS line after LOCUS $locus"
                if defined $locus;
            $locus = $1;
            print STDERR "Found locus $locus\n";
        } elsif ($line eq "ORIGIN") {
            die "Unexpected ORIGIN line"
                if !defined $locus;
            die "Duplicate ORIGIN line in $locus" if defined $seqsofar;
            $seqsofar = "";
        } elsif (defined $seqsofar && $line =~ m/^ *[0-9]+ ([A-Za-z ]+)$/) {
            $seqsofar .= $1;
        } elsif ($line eq "//") {
            die "Ending section with no locus" unless defined $locus;
            die "Ending locus $locus with no ORIGIN" unless defined $seqsofar;
            die "No sequence for locous $locus" if $seqsofar eq "";
            $seqsofar =~ s/ //g;
            $seqsofar = uc($seqsofar);
            my @seqs = $seqsofar =~ /.{1,60}/g;
            print ">$locus", "\n";
            foreach my $seq (@seqs) { print $seq, "\n"; }
            $locus = undef;
            $seqsofar = undef;
        }
    }
}
