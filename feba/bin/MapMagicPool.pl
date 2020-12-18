#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable() and reverseComplement()
use FindGene; # for LocationToGene()

my $gdir = "g";
my $minconf = 50;

my $usage = <<END
MapMagicPool.pl -org nickname -library library-name
                -magic magic_construct_barcodes
                file1.mapping ... fileN.mapping

Given mapping files (from MapTnSeq.pl or RunTnSeqLocal.pl) and a table
of (BarSeq) barcode to construct, analyzes the behavior of well
represented constructs.

Writes to g/nickname/library_magic.* -- a hiconf file with the
high-confidence mappings, and a bias file with statistics for the
popular constructs. Both are tab-delimited.

magic_construct_barcodes should be tab-delimited with the fields
construct and barcode (in the BarSeq orientation).

In the hiconf file, barcode is the sequence in the TnSeq direction and
rcbarcode is the sequence in the BarSeq direction. locusId is blank if
the insertion is not in a gene.

The bias file has the fields:
name -- the construct, as in the input table
nReads -- total reads for the constructs' barcodes (either
          past-the-end of the transposon or mapped into genome)
nPstEnd -- past-the-end reads for this construct
fPstEnd -- fraction of reads for intact vector
nLoc -- number of high-confidence insertions
nInGene -- number of insertions within genes
fStrAgr -- fraction of time gene strand = insertion strand
        (for a good vector, will be very close to 50%)
nGenesH -- number of different genes with insertions
The remaining statistics are computed while considering only
genes that are hit (with at least one insertion):
meanLoc -- average of #insertions per gene
medLoc -- median of #insertions per gene
locBias -- mean/median #insertions per gene
meanRd -- average of #reads per gene
medRd -- median of #reads per gene
rdBias -- mean/median #reads per gene

Optional arguments:
   -gdir $gdir -- parent of the nickname/ directory
   -minconf $minconf -- minimum number of high-confidence locations
	for a construct before it is reported.
END
    ;

sub mean;
sub median;

{
    my ($org,$libname,$magicfile);
    die $usage unless
        GetOptions('org=s' => \$org,
                   'library=s' => \$libname,
                   'magic=s' => \$magicfile,
                   'gdir=s' => \$gdir,
                   'minconf=i' => \$minconf )
        && defined $org
        && defined $libname
        && defined $magicfile;
    die "No such directory: $gdir\n" unless -d $gdir;
    die "No such directory: $gdir/$org\n" unless -d "$gdir/$org";
    my $outpre = "$gdir/$org/${libname}_magic";

    my $genesfile = "$gdir/$org/genes.GC";
    die "No such file: $genesfile\n" unless -e $genesfile;
    my @genes = ReadTable($genesfile, [ 'locusId', 'scaffoldId', 'strand', 'begin', 'end' ]);
    my %geneScaffolds = map { $_->{scaffoldId} => 1 } @genes;
    my %genesSorted = (); # scaffold to list of genes sorted by begin
    foreach my $gene (@genes) {
	push @{ $genesSorted{$gene->{scaffoldId}} }, $gene;
    }
    foreach my $scaffold (keys %genesSorted) {
	my @sorted = sort { $a->{begin} <=> $b->{begin} } @{ $genesSorted{$scaffold} };
	$genesSorted{$scaffold} = \@sorted;
    }
    CheckGeneLocations(\%genesSorted); # writes to STDERR
    my %locusStrand = map { $_->{locusId} => $_->{strand} } @genes;

    die "No such file: $magicfile\n" unless -e $magicfile;
    my @mappingfiles = @ARGV;
    die "No mapping files specified\n" unless scalar(@mappingfiles) > 0;
    foreach my $file (@mappingfiles) {
        die "No such file: $file\n" unless -e $file;
    }

    my @constructs = ReadTable($magicfile, ['barcode', 'construct']);
    print STDERR "Read " . scalar(@constructs) . " constructs from $magicfile\n";
    my %rcToConstruct = map {
        reverseComplement( $_->{barcode} ) => $_->{construct}
    } @constructs;
    die "Some barcodes are reused in $magicfile\n"
        unless scalar(keys %rcToConstruct) == scalar(@constructs);

    my %constructPastEnd = ();
    my %constructMapped = ();
    my %count = (); # barcode => join("\t", scaffold, pos, strand) => count
    my $nReads = 0;
    my $nMatch = 0; # reads for known barcodes
    foreach my $file (@mappingfiles) {
        open(FILE, "<", $file) || die "Cannot read $file";
        while(my $line = <FILE>) {
            chomp $line;
            $nReads++;
            my ($read,$barcode,$scaffold,$pos,$strand,$uniq,$sbeg,$send,$bits,$identity) = split /\t/, $line, -1;
            die "Invalid line\n$line\nin $file"
                unless defined $identity;
            next unless exists $rcToConstruct{$barcode};
            my $con = $rcToConstruct{$barcode};
            if ($scaffold eq "pastEnd") {
                $constructPastEnd{$con}++;
                $nMatch++;
            } else {
                die "Invalid line\n$line\nin $file"
                    unless $pos =~ m/^\d+$/
                    && ($strand eq "+" || $strand eq "-")
                    && $scaffold ne "";
                $count{$barcode}{ join("\t",$scaffold,$pos,$strand) }++ if $uniq eq "1";
                $constructMapped{$con}++;
                $nMatch++;
            }
        }
        close(FILE) || die "Error reading $file";
    }
    die "No reads\n" if $nReads == 0;
    print STDERR sprintf("Reads %d Matching known barcodes %d (%.3f%%)\n",
                         $nReads, $nMatch, (100.0*$nMatch)/$nReads);

    # construct => list of [ barcode, scaffold, pos, strand, nreads, locusId or "", locusStrand or "" ]
    # nreads > 1 only
   open(HICONF, ">", "$outpre.hiconf")
       || die "Cannot write to $outpre.hiconf";
   print HICONF join("\t", "construct", "barcode", "rcbarcode",
                     "scaffold", "pos", "strand", "n", "locusId", "locusStrand")."\n";

    my %hiconf = ();
    while (my ($barcode, $hash) = each %count) {
        while (my ($poscode, $n) = each %$hash) {
            next unless $n > 1;
            my ($scaffold,$pos,$strand) = split /\t/, $poscode;
            die if !defined $strand;
            my $con = $rcToConstruct{$barcode};
            die unless defined $con;
	    my ($locusId, $f) = &LocationToGene($scaffold, $pos, \%genesSorted);
            my $locusStrand = $locusStrand{$locusId} || "";
            $locusId = "" if !defined $locusId;
            push @{ $hiconf{$con} }, [ $barcode, $scaffold, $pos, $strand, $n, $locusId, $locusStrand ];
            print HICONF join("\t", $con, $barcode, reverseComplement($barcode),
                              $scaffold, $pos, $strand, $n, $locusId, $locusStrand)."\n";
        }
    }
    close(HICONF) || die "Error writing to $outpre.hiconf";
    print STDERR "Wrote $outpre.hiconf\n";

    # And report bias of each construct with enough high-confidence locations
    open(BIAS, ">", "$outpre.bias") || die "Cannot write to $outpre.bias";
    print BIAS join("\t", qw{name nReads nPstEnd fPstEnd nLoc nInGene fStrAgr nGenesH meanLoc medLoc locBias meanRd medRd rdBias})."\n";
    foreach my $con (sort keys %hiconf) {
        my $list = $hiconf{$con};
        next unless scalar(@$list) >= $minconf && scalar(@$list) > 0;
        my $nInGene = 0;
        my $nStrandAgree = 0;
        my %geneNStrains = ();
        my %geneNReads = ();
        foreach my $row (@$list) {
            my ($barcode,$scaffold,$pos,$strand,$n, $locusId, $locusStrand) = @$row;
            die unless $n > 1;
            next unless $locusId ne "";
            $nInGene++;
            die unless defined $locusStrand && ($locusStrand eq "+" || $locusStrand eq "-");
            $nStrandAgree++ if $strand eq $locusStrand;
            $geneNStrains{$locusId}++;
            $geneNReads{$locusId} += $n;
        }
        my $pastEnd = $constructPastEnd{$con} || 0;
        my @stats = ();
        if ($nInGene == 0) {
            @stats = ("","","","","","");
        } else {
            my $meanLoc = mean(values %geneNStrains);
            my $medLoc = median(values %geneNStrains);
            my $meanReads = mean(values %geneNReads);
            my $medReads = median(values %geneNReads);
            @stats = ( sprintf("%.1f", $meanLoc), sprintf("%.1f",$medLoc),
                       sprintf("%.2f", $meanLoc/$medLoc),
                       sprintf("%.1f", $meanReads), sprintf("%.1f", $medReads),
                       sprintf("%.2f", $meanReads/$medReads) );
        }
        print BIAS join("\t",
                        $con,
                        $constructMapped{$con} + $pastEnd, $pastEnd,
                        sprintf("%.4f",$pastEnd/($constructMapped{$con} + $pastEnd)),
                        scalar(@$list),
                        $nInGene,
                        sprintf("%.4f",$nStrandAgree/$nInGene),
                        scalar(keys %geneNStrains),
                        @stats
            )."\n";
    }
    close(BIAS) || die "Error writing $outpre.bias";
    print STDERR "Wrote $outpre.bias\n";
}

sub mean {
    my (@in) = @_;
    @in = grep { defined $_ } @in;
    return undef unless scalar(@in) > 0;
    my $sum = 0;
    foreach my $val (@in) {
        $sum += $val;
    }
    return $sum / scalar(@in);
}

sub median {
    my (@in) = @_;
    @in = grep { defined $_ } @in;
    return undef unless scalar(@in) > 0;
    @in = sort { $a <=> $b } @in;
    my $midI = int(scalar(@in)/2);
    if (scalar(@in) % 2 == 0) {
        # i.e. for 4, use indexes 1 and 2, or the 2nd and 3rd elements
        return ($in[$midI-1] + $in[$midI])/2;
    }
    #else
    return $in[$midI];
}
