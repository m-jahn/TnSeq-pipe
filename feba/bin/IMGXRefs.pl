#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()

my $gdir = "g";
my $usage = <<END
Usage: IMGXRefs.pl [ -g $gdir ] nickname1 ... nicknameN > xrefs.tab
    g/nickname/IMG.tab files should contain the gene table exported from
    JGI IMG, with fields named "gene_oid" and "Locus Tag". If the values in
    "Locus Tag" match to sysName in the genes table (in g/nickname/genes.tab),
    then entries are written to the output. Only proteins (genes.tab type=1) are considered.
    If the IMG.tab file is missing, or if none of the locus tags match, then a warning
    is issued to standard error.
END
    ;

{
    (GetOptions('g=s' => \$gdir) && @ARGV > 0) || die $usage;
    my @orgs = @ARGV;
    my @orgs2 = (); # orgs to try to match for

    foreach my $org (@orgs) {
        die "No such directory: $gdir/$org" unless -d "$gdir/$org";
        die "No genes table (genes.tab) for $org" unless -e "$gdir/$org/genes.tab";
        if (-e "$gdir/$org/IMG.tab") {
            push @orgs2, $org;
        } else {
            print STDERR "No IMG.tab file for $org\n";
        }
    }
    my $nSuccess = 0;
    print join("\t", qw{orgId locusId xrefDb xrefId})."\n";
    foreach my $org (@orgs2) {
        my @genes = &ReadTable("$gdir/$org/genes.tab",qw{locusId sysName type});
        my %sysNameToLocusId = ();
        foreach my $gene (@genes) {
            $sysNameToLocusId{$gene->{sysName}} = $gene->{locusId} if $gene->{type} eq "1";
        }
        my @img = &ReadTable("$gdir/$org/IMG.tab",("gene_oid","Locus Tag"));
        my $nMatch = 0;
        my $lastid = "";
        foreach my $row (@img) {
            my $locustag = $row->{"Locus Tag"};
            next if $locustag eq $lastid;
            $lastid = $locustag;
            my $locusId = $sysNameToLocusId{$locustag};
            if ($locustag ne "" && defined $locusId && $row->{gene_oid} =~ m/^\d+$/) {
                die "$org $locustag empty locusId" if $locusId eq "";
                print join("\t", $org, $locusId, "IMG", $row->{gene_oid})."\n";
                $nMatch++;
            }
        }
        if ($nMatch > 0) {
            $nSuccess++;
        } else {
            print STDERR "No matches by locus tag for $org\n";
        }
    }
    print STDERR "Found matches for $nSuccess of " . scalar(@orgs2) . " genomes with IMG.tab files\n";
}
