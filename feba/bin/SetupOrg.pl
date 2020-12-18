#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable

my $usage = <<END
SetupOrg.pl [ -prefix ABC_ ] -gbk file.gbk [ -name NickName | -out outdir ]
	(where file.gbk is from IMG, RAST, or RefSeq)
SetupOrg.pl -name NickName -gff gff_file -fna genome.fna [ -faa prot.faa ] [ -name NickName | -out outdir ]
	(for assemblies from the JGI)

Creates the output directory (if it does not exist already) and within
it, it creates the genome sequence file (genome.fna), the genes files
(genes.tab, genes.GC), and the protein sequence file (aaseq).

The default output directory is g/NickName

If using a genbank file, please make sure that it includes all scaffolds.

Use the -aaseq option to write the aaseq file only

Use the -prefix argument to create locusIds like ABC_ not GFFnnnn
(this is mostly relevant if reading genbank files from RAST).

Dependencies: Many of the scripts that this calls require BioPerl.
Side effects: deletes aaseq2, pfam.tab, and tigrfam.tab from the output directory (if they exist)
END
    ;

# Given the genes table and the input faa file, make sure that the entries
# in the faa file correspond to known genes, and change the type as necessary
# to match presence/absence in the faa file.
#
# Also, cleans up the faa file to have locusId as the id and writes it to the 3rd argument
sub CheckAASeq($$$);

{
    my ($name, $gbkFile, $aaseqOnly, $fnaFile, $gffFile, $faaFile, $prefix, $outdir);
    $prefix = "GFF";

    die $usage unless GetOptions('aaseq' => \$aaseqOnly,
				 'name=s' => \$name,
				 'gbkFile=s' => \$gbkFile,
				 'gffFile=s' => \$gffFile,
				 'fnaFile=s' => \$fnaFile,
                                 'faa=s' => \$faaFile,
				 'prefix=s' => \$prefix,
				 'out=s' => \$outdir )
	&& @ARGV==0;
    die $usage unless (defined $gbkFile) xor (defined $gffFile);
    die $usage unless (defined $gffFile && defined $fnaFile) || (!defined $gffFile && !defined $fnaFile);
    die $usage unless (defined $name) xor (defined $outdir);
    die "Cannot specify -faa unless using -gff as well\n"
        if defined $faaFile && !defined $gffFile;
    die "No such file: $faaFile"
        if defined $faaFile && ! -e $faaFile;

    die "Invalid prefix $prefix" unless $prefix =~ m/^[A-Za-z][A-Za-z0-9_]*$/;

    if (!defined $outdir) {
	die "No g directory" unless -d "g";
	$outdir = "g/$name";
    }
    if (! -d $outdir) {
	mkdir($outdir) || die "Cannot mkdir $outdir";
	system("chmod","g+w",$outdir);
    }
    if (defined $gbkFile) {
	die "No such file: $gbkFile" unless -e $gbkFile;
	print STDERR "Creating $outdir from $gbkFile\n";
    }
    print STDERR "Creating aaseq file only\n" if defined $aaseqOnly;

    # these files should be created from aaseq, so delete them before writing to aaseq
    unlink("$outdir/aaseq2");
    unlink("$outdir/pfam.tab");
    unlink("$outdir/tigrfam.tab");

    if (defined $gbkFile && !defined $aaseqOnly) {
        die "Please install genbank2gff.pl in $Bin and make it executable\n"
            . "This script is available from Ian Holmes' gfftools repository\n"
            . "See https://github.com/ihh/gfftools\n"
            unless -x "$Bin/genbank2gff.pl";
	system("$Bin/gbkToSeq.pl $gbkFile > $outdir/genome.fna") == 0
            || system("$Bin/gbkToSeq2.pl $gbkFile > $outdir/genome.fna") == 0
            || die "Both gbkToSeq.pl and gbkToSeq2.pl failed";
	system("$Bin/genbank2gff.pl $gbkFile > $outdir/genes.gff") == 0 || die "genbank2gff.pl failed";
	system("$Bin/gffToGenes.pl -prefix $prefix < $outdir/genes.gff > $outdir/genes.tab") == 0 || die "gffToGenes.pl failed";
    }
    if (defined $gbkFile) {
	my $code = system("$Bin/gbkToFaa.pl $outdir/genes.tab $gbkFile > $outdir/aaseq");
	if ($code != 0) {
	    print STDERR "Warning: gbkToFaa.pl failed (common with RAST genbank files) -- translating the nt sequences instead\n";
	    system("$Bin/genesTabTranslation.pl $outdir/genes.tab $outdir/genome.fna > $outdir/aaseq") == 0
		|| die "genesTabTranslation.pl failed";
	}
    }

    if (defined $gffFile && !defined $aaseqOnly) {
	system("cp", $fnaFile, "$outdir/genome.fna") == 0 || die $!;
	system("$Bin/gffToGenes.pl -prefix $prefix < $gffFile > $outdir/genes.tab") == 0 || die "gffToGenes.pl failed";
    }
    if (defined $gffFile) {
        if (defined $faaFile) {
            CheckAASeq("$outdir/genes.tab", $faaFile, "$outdir/aaseq");
        } else {
            system("$Bin/genesTabTranslation.pl $outdir/genes.tab $outdir/genome.fna > $outdir/aaseq") == 0
                || die "genesTabTranslation.pl failed";
        }
    }

    if (!defined $aaseqOnly) {
	system("$Bin/RegionGC.pl $outdir/genome.fna $outdir/genes.tab > $outdir/genes.GC") == 0 || die "RegionGC.pl failed";
    }
}

sub CheckAASeq($$$) {
    my ($geneTabFile, $faaIn, $faaOut) = @_;
    my @req = qw{locusId sysName type};
    my @geneCol = ReadColumnNames($geneTabFile);
    my @genes = ReadTable($geneTabFile, \@req);
    die "No genes in $geneTabFile" unless @genes > 0;
    my %idToGene = map { $_->{locusId} => $_ } @genes;
    my %tagToGene = map { $_->{sysName} => $_ } @genes;

    my %seqs = ();
    open(IN, "<", $faaIn) || die "Cannot read $faaIn";
    my $cur = undef;
    while(<IN>) {
        if (m/^>(\S+) (\S+)/) {
            my ($name1,$name2) = ($1,$2);
            if (exists $idToGene{$name1}) {
                $cur = $name1;
            } elsif (exists $idToGene{$name2}) {
                $cur = $name2;
            } elsif (exists $tagToGene{$name1}) {
                $cur = $tagToGene{$name1}{locusId};
            } elsif (exists $tagToGene{$name2}) {
                $cur = $tagToGene{$name2}{locusId};
            } else {
                die "Unrecognized identifier in $faaIn: $name1 $name2\n";
            }
            die "Duplicate identifier $cur (from $name1 $name2) in $faaIn\n"
                if exists $seqs{$cur};
            $seqs{$cur} = "";
        } elsif (m/^>(\S+)/) {
            my $name = $1;
            if (exists $idToGene{$name}) {
                $cur = $name;
            } elsif (exists $tagToGene{$name}) {
                $cur = $tagToGene{$name}{locusId};
            } else {
                die "Unrecognized identifier in $faaIn: $name\n";
            }
            die "Duplicate identifier $cur (from $name) in $faaIn\n"
                if exists $seqs{$cur};
            $seqs{$cur} = "";
        } elsif (m/^>/) {
            die "Cannot parse fasta header line\n$_";
        } else {
            chomp;
            die "No header line at start of fasta file $faaIn" unless defined $cur;
            s/\s//g; # remove whitespace
            s/[*]/X/g; # sometimes used as an ambiguous or stop character
            die "Invalid characters in protein fasta line $_"
                unless m/^[A-Z]*$/; # allow empty line!
            $seqs{$cur} .= $_;
        }
    }
    close(IN) || die "Error reading $faaIn";
    die "No sequences in $faaIn" unless scalar(keys %seqs) > 0;
    print STDERR "Read sequences for " . scalar(keys %seqs) . " proteins from $faaIn\n";
    while (my ($locusId,$aaseq) = each %seqs) {
        die "Empty sequence for $locusId in $faaIn" if $seqs{$locusId} eq "";
    }

    foreach my $gene (@genes) {
        my $locusId = $gene->{locusId};
        my $sysName = $gene->{sysName};
        if ($gene->{type} eq "1" && !exists $seqs{$locusId}) {
            print STDERR "No sequence for protein $locusId $sysName -- changing to pseudogene\n";
            $gene->{type} = "7";
        } elsif ($gene->{type} ne "1" && exists $seqs{$locusId}) {
            print STDERR "Found protein sequence for $locusId $sysName -- changing from type $gene->{type} to 1 (protein)\n";
            $gene->{type} = "1";
        }
    }

    # save the protein sequences
    open(FAA, ">", $faaOut) || die "Cannot write to $faaOut";
    foreach my $id (sort keys %seqs) {
        my $seq = $seqs{$id};
        print FAA ">$id\n$seq\n";
    }
    close(FAA) || die "Error writing to $faaOut";
    print STDERR "Wrote $faaOut\n";

    # save the revised gene table
    open(GENES, ">", $geneTabFile) || die "Cannot write to $geneTabFile";
    print GENES join("\t", @geneCol)."\n";
    foreach my $gene (@genes) {
        my @out = ();
        foreach my $col (@geneCol) {
            die "Missing field $col in gene" unless exists $gene->{$col};
            push @out, $gene->{$col};
        }
        print GENES join("\t", @out)."\n";
    }
    close(GENES) || die "Error writing to $geneTabFile";
    print STDERR "Wrote corrected $geneTabFile\n";
}
