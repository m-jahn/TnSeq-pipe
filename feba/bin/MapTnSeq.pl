#!/usr/bin/perl -w
# Given TnSeq data from a library of transposon insertions with random barcodes, for each usable read,
# identify the barcode and the location in the genome.
#
# -- Morgan Price, Arkin group, Lawrence Berkeley National Lab

use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $minQuality = 10; # every nucleotide in a barcode must be at least this quality
my $flanking = 5; # number of nucleotides on each side that must match
my $wobbleAllowed = 2; # uncertainty in location of barcode or end of transposon, on either side of expectation
my $tmpdir = defined $ENV{TMPDIR} ? $ENV{TMPDIR} : "/tmp";
my $minIdentity = 90; # minimum %identity for mapping to genome or past-end
my $minScore = 15; # minimum score for mapping to genome or past-end
my $delta = 5; # minimum difference in score for considering a hit unique
my $debug = undef;
# Parameters for BLAT. Thanks to Judy Savitskaya for these additions
my $tileSize = 11; # size of an alignment tile
my $stepSize = 11; # distance between the starting bases of alignment tiles (will overlap if stepSize<tileSize)

# Given BLAT rows (as list of lists) and hitsPastEnd hash (as reference),
# output the mapping for the read, or not (if "trumped" by hit-past-end).
# Updates hitsPastEnd for the read to score=0 if genomic hit trumps hit-past-end.
sub HandleGenomeBLAT($$);

# Global so that HandleGenomeBLAT() can use them
my $nMapped = 0;
my $nMapUnique = 0;
my $nPastEndIgnored = 0; # weak hit to past-end ignored
my $nPastEndTrumps = 0; # hit to past-end (close to) as good as hit to genome
my %nameToBarcode = ();

# Note -- should also add an option for tweaking tileSize and stepSize? 12/4 seems better than the default 11/11?

my $usage = <<END
Usage: MapTnSeq.pl [ -debug ] [ -limit maxReads ] [ -minQuality $minQuality ]
            [-flanking $flanking]  [ -wobble $wobbleAllowed ]
            [ -minIdentity $minIdentity ] [ -minScore $minScore ] [ -delta $delta ]
            [ -tileSize $tileSize ] [ -stepSize $stepSize ]
            [-tmpdir $tmpdir ] [ -unmapped saveTo ] [ -trunc saveTo ]
            -genome fasta_file -model model_file -first fastq_file > output_file

    The fastq file should have phred+33 ("sanger") encoding of quality scores
    (as in MiSeq or 2012+ HiSeq). If it is named .gz, it will be gunzipped before reading.
    If it contains paired-end reads, the second read (name matching " 2:") will be ignored.

    The model file contains 1-2 lines -- the first shows what a
    typical read should look like up until the junction with the genome, e.g.
    
    nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT

    All characters in the model read must be ACGT except for the
    optional block of n\'s at the front and a block of Ns which
    represents the barcode.

    The second line in the model file is optional and contains the
    sequence "past the end" of the transposon that might arise from
    residual intact plasmid.

    This script does not handle sample multiplexing.

    The output file is tab-delimited and contains, for each usable
    read, the read name, the barcode, which scaffold the insertion
    lies in, the position of the insertion, the strand that the read
    matched, a boolean flag for if this mapping location is unique or
    not, the beginning and end of the hit to the genome in the read
    after trimming the transposon sequence, the bit score, and the
    %identity.

    minQuality specifies the minimum quality score for each character
    in the barcode.  flanking specifies the minimum number of
    nucleotides on each side that must match exactly.

    tileSize and stepSize are parameters for BLAT.  If the portion of
    the read after the junction is short, try a lower stepSize.  Also
    see BLAT's documentation.

    Use -unmapped or -trunc to save the portion of the reads past the
    junction. Only reads with barcodes and junctions are
    saved. -unmapped writes only the unmapped reads, in fasta format,
    with the barcode and the read name as the definition line. -trunc
    writes the remaining part of all of the reads with junctions, in
    fastq format, and it appends :barcode to the end of the read name
    (but before the space).

END
    ;

sub FindBarcode($$$$$);
sub FindModelEnd($$$);
sub FindSubstr($$$$);
sub BLAT8($$$$$$$$); # BLAT to a blast8 format file

{
    my $limit = undef;
    my $fastqFile = undef;
    my $genomeFile = undef;
    my $modelFile = undef;
    my $blatcmd = -e "$Bin/blat" ? "$Bin/blat" : "blat";
    my $unmappedFile = undef;
    my $truncFile = undef;
    
    my $minGenomeId = 90;

    (GetOptions('debug' => \$debug, 'limit=i' => \$limit, 'minQuality=i' => \$minQuality,
		'unmapped=s' => \$unmappedFile,
                'trunc=s' => \$truncFile,
                'flanking=i' => \$flanking, 'minIdentity=i' => \$minIdentity, 'minScore=i' => \$minScore,
                'delta=i' => \$delta,
                'tileSize=i' => \$tileSize,'stepSize=i' => \$stepSize,
                'tmpdir=s' => \$tmpdir,
                'blat=s' => \$blatcmd,
                'wobble=s' => \$wobbleAllowed,
                'genome=s' => \$genomeFile, 'model=s' => \$modelFile, 'first=s' => \$fastqFile)
     && @ARGV==0) || die $usage;
    die $usage unless defined $modelFile && defined $fastqFile;

    die "Cannot read $genomeFile" unless -r $genomeFile;
    die "Cannot read $fastqFile" unless -r $fastqFile;
    die "Not a directory: $tmpdir" unless -d $tmpdir;
    die "minScore must be at least 10" if $minScore < 10;
    die "tileSize must be at least 7" if $tileSize < 7;
    die "stepSize must be at least 1" if $stepSize < 1;
    die "Invalid minimum identity" if $minIdentity < 50 || $minIdentity > 100;
    die "delta cannot be negative\n" if $delta < 0;

    open(MODEL, "<", $modelFile) || die "Cannot read $modelFile";
    my $model = <MODEL>;
    $model =~ s/[\r\n]+$//;
    die "Invalid model: $model" unless $model =~ m/^n*[ACGT]+N+[ACGT]+$/;
    my $pastEnd = <MODEL>;
    if (defined $pastEnd) {
	chomp $pastEnd;
	die "Invalid past-end sequence: $pastEnd" unless $pastEnd =~ m/^[ACGT]+$/;
    }	
    close(MODEL) || die "Error reading $modelFile";

    my $barcodeStart = index($model, "N");
    my $barcodeEnd = rindex($model, "N");
    my $barcodeLen = $barcodeEnd-$barcodeStart+1;
    print STDERR "Parsed model $modelFile\n";
    print STDERR "Barcodes of length $barcodeLen, expected transposon region of " . length($model);
    print STDERR ", no pastEnd" if !defined $pastEnd;
    print STDERR "\n";

    my $pipe = 0;
    if ($fastqFile =~ m/[.]gz$/) {
        $pipe = 1;
        open(FASTQ, '-|', 'zcat', $fastqFile) || die "Cannot run zcat on $fastqFile";
        print STDERR "Reads (gzipped) from $fastqFile\n";
    } elsif ($fastqFile =~ m/[.zip]$/) {
        $pipe = 1;
        open(FASTQ, "7za -so e $fastqFile | zcat |") || die "Cannot run 7za and zcat on $fastqFile";
        print STDERR "Reads (zipped) from $fastqFile\n";
    } else {
        open(FASTQ, "<", $fastqFile) || die "Cannot read $fastqFile";
        print STDERR "Reads (raw) from $fastqFile\n";
    }

    my %nameToHits = (); # list of scaffold, position, strand, match score
    my $nLong = 0;

    my $rand = rand();
    my $tmpFna = $tmpdir . "/MapTnSeq_" . $$ . "_$rand.fna";
    open(TMPFNA, ">", $tmpFna) || die "Cannot write to $tmpFna";

    my $trunc_fh;
    if ($truncFile) {
      open($trunc_fh, ">", $truncFile) || die "Cannot write to $truncFile";
    }

    my $nReads = 0;
    my $nTryToMap = 0;
    # read name => 1 if mapped or pastEnd, 0 otherwise
    # to save memory, only fill out with the -unmapped option
    my %mapnames  = ();

    # Find barcodes and end of transposon and write remaining sequence to TMPFNA
    while(!defined $limit || $nReads < $limit) {
        my $name = <FASTQ>;
        (last, next) unless $name;
        chomp $name;
        die "Sequence name line does not start with @" unless $name =~ m/^@/;

        my $seq = <FASTQ>;
        chomp $seq;
        die "Sequence line is empty"
            unless defined $seq && length($seq) > 0;
        die "Sequence line contains invalid characters: $seq"
            unless $seq =~ m/^[A-Z]+$/;
        my $line = <FASTQ>;
        die "Third line does not start with +"
            unless defined $line && $line =~ m/^[+]/;
        my $quality = <FASTQ>;
        chomp $quality;
        die "Quality line is of wrong length"
            unless defined $quality && length($quality) == length($seq);

        next if $name =~ m/^\S+ 2:/; # ignore second side of paired-end reads
        $nReads++;

        # short sequences are unmappable
        next unless length($seq) >= length($model) + $minScore;
        $nLong++;

        my ($barcode,$obsStart) = FindBarcode($seq,$quality,$model,$barcodeStart,$barcodeEnd);
        next unless defined $barcode;

        my $shortname = $name; $shortname =~ s/ .*$//;
        die "Duplicate read name: $shortname" if exists $nameToBarcode{$shortname};
        $nameToBarcode{$shortname} = $barcode;

        my $transposonEnd = FindModelEnd($seq,$model,$obsStart - $barcodeStart);
        if (defined $transposonEnd && length($seq) >= $transposonEnd + $minScore) {
            print STDERR "Try to map $shortname\n" if $debug;
            my $inGenome = substr($seq, $transposonEnd+1);
            print TMPFNA ">$shortname\n$inGenome\n";
            $nTryToMap++;
	    $mapnames{$shortname} = 0 if defined $unmappedFile;
            if ($truncFile) {
              my @words = split / /, $name;
              $words[0] .= ":$barcode";
              print $trunc_fh join(" ", @words), "\n",
                substr($seq, $transposonEnd), "\n",
                "+", "\n",
                substr($quality, $transposonEnd), "\n";
            }
        }
    }

    # Prematurely closing the pipe gives an error
    close(FASTQ) || ($pipe && defined $limit) || die "Error reading from $fastqFile: $!";
    close(TMPFNA) || die "Error writing to $tmpFna";
    print STDERR "Read $nReads reads\n";
    if ($truncFile) {
      close($trunc_fh) || die "Error writing to $truncFile";
    }

    if ($nTryToMap == 0) {
      print STDERR "None of the reads are candidates for mapping (none match the model and are long enough)\n";
      unlink($tmpFna);
      exit(0);
    }

    my %hitsPastEnd = (); # read to score of match to past-end sequence

    if (defined $pastEnd) {
	# Map to past-end-of transposon
	my $endFna = $tmpdir . "/MapTnSeq_end" . $$ . "_" . "_$rand.fna";
	open(END, ">", $endFna) || die "Cannot write to $endFna";
	print END ">pastend\n$pastEnd\n";
	close(END) || die "Error writing to $endFna";
	my $blat8 = BLAT8($tmpFna, $endFna, $tmpdir, $blatcmd, $minScore, $minIdentity,$tileSize,$stepSize);
	print STDERR "Parsing past-end hits to $blat8\n" if defined $debug;
	open(BLAT, "<", $blat8) || die "Cannot read $blat8";
	while(<BLAT>) {
	    chomp;
	    my @F = split /\t/, $_;
	    my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;
	    $hitsPastEnd{$query} = $score unless exists $hitsPastEnd{$query} && $hitsPastEnd{$query} > $score;
	    $mapnames{$query} = 1 if defined $unmappedFile;
	}
	close(BLAT) || die "Error reading $blat8";
	unlink($blat8) unless defined $debug;
    }

    # Map to the genome
    my $blat8 = BLAT8($tmpFna, $genomeFile, $tmpdir, $blatcmd, $minScore, $minIdentity, $tileSize, $stepSize);
    print STDERR "Parsing $blat8\n" if defined $debug;
    open(BLAT, "<", $blat8) || die "Cannot read $blat8";
    my @lines = ();
    while(<BLAT>) {
        chomp;
        my @F = split /\t/, $_;
        my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @F;
	$mapnames{$query} = 1 if defined $unmappedFile;
        if (@lines == 0 || $query eq $lines[0][0]) {
            push @lines, \@F;
        } else {
            HandleGenomeBLAT(\@lines, \%hitsPastEnd);
            @lines = \@F;
        }
    }
    HandleGenomeBLAT(\@lines, \%hitsPastEnd);

    close(BLAT) || die "Error reading $blat8";
    unlink($blat8) unless defined $debug;

    if (defined $unmappedFile) {
	my $nUnmapped = 0;
	open(UNMAPPED, ">", $unmappedFile) || die "Cannot write to $unmappedFile";
	open(FNA, "<", $tmpFna) || die "Cannot read $tmpFna";
	while(my $header = <FNA>) {
	    chomp $header;
	    $header =~ m/^>(.*)$/ || die "Cannot parse $header";
	    my $name = $1;
	    my $seq = <FNA>;
	    chomp $seq;
	    $seq =~ m/^[ACGTN]+$/ || die "Cannot parse $seq in $tmpFna";
	    if ($mapnames{$name} == 0) {
		print UNMAPPED ">$nameToBarcode{$name} $name\n$seq\n";
		$nUnmapped++;
	    }
	}
	close(FNA) || die "Error reading $tmpFna";
	close(UNMAPPED) || die "Error writing to $unmappedFile";
	print STDERR "Wrote $nUnmapped unmapped reads to $unmappedFile in fasta format\n";
    }

    unlink($tmpFna) unless defined $debug;

    # print out hits-past-end
    while (my ($read,$score) = each %hitsPastEnd) {
        next unless $score > 0; # score=0 means it was ignored above
        print join("\t", $read, $nameToBarcode{$read}, "pastEnd", "", "", "", "", "", "", "")."\n";
    }

    print STDERR "Reads processed $nReads Long-enough $nLong Barcodes found " . scalar(keys %nameToBarcode) . "\n";
    print STDERR "Mapping attempted for $nTryToMap Mapped $nMapped Uniquely $nMapUnique\n";
    print STDERR "Hits past end of transposon: " . (scalar(keys %hitsPastEnd) - $nPastEndIgnored) .
        " plus $nPastEndIgnored weak/ambiguous; trumped hit to genome $nPastEndTrumps times\n";
    print STDERR sprintf("Proportions: Long-enough %.3f Barcode %.3f Attempted %.3f Mapped %.3f Past-end %.3f\n",
			 $nLong/$nReads,
			 scalar(keys %nameToBarcode)/$nReads,
			 $nTryToMap/$nReads,
			 $nMapped/$nReads,
			 (scalar(keys %hitsPastEnd) - $nPastEndIgnored)/$nReads)
	if $nReads > 0;
}

sub FindSubstr($$$$) {
    my ($subseq, $seq, $expAt, $wobble) = @_;
    my $len = length($seq);
    for (my $i = $expAt - $wobble; $i <= $expAt+$wobble; $i++) {
        return($i) if $i >= 0 && $i < $len && substr($seq, $i, length($subseq)) eq $subseq;
    }
    return undef;
}

sub FindBarcode($$$$$) {
    my ($seq,$quality,$model,$expStart,$expEnd) = @_;
    my $pre = substr($model, $expStart - $flanking, $flanking);
    my $post = substr($model, $expEnd + 1, $flanking);
    my $preLoc = FindSubstr($pre, $seq, $expStart - $flanking, $wobbleAllowed);
    return undef unless defined $preLoc;
    my $postLoc = FindSubstr($post, $seq, $expEnd+1, $wobbleAllowed);
    return undef unless defined $postLoc;

    # positions of 1st and last nucleotides of barcode
    my $start = $preLoc + $flanking;
    my $end = $postLoc-1;

    unless ($end-$start eq $expEnd-$expStart) {
        print STDERR "Wrong barcode length: $start $end not $expStart $expEnd in $seq\n" if defined $debug;
        return undef;
    }
    my $barcode = substr($seq, $start, $end-$start+1);

    if ($minQuality > 0) {
        # the sanger code for !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
        my $barqual = substr($quality, $start, $end-$start+1);
        # quality score is encoded as 33+
        my @scores = map { $_ - 33 } unpack("%C"x20, $barqual);
        foreach my $score (@scores) {
            die "Invalid score $score from barcode $barcode quality $barqual" if $score < 0 || $score > 100;
            if ($score < $minQuality) {
                print STDERR "Low quality $score for barcode $barcode in $seq\n" if defined $debug;
                return undef;
            }
        }
    }
    return($barcode,$start);
}

sub FindModelEnd($$$) {
    my ($seq,$model,$expOff) = @_;
    my $expEnd = length($model) + $expOff;
    my $at = FindSubstr(substr($model,length($model)-$flanking,$flanking), $seq, $expEnd-$flanking, $wobbleAllowed);
    if (!defined $at) {
        print STDERR "No end sequence " . substr($model,length($model)-$flanking,$flanking) . " near position " . ($expEnd-$flanking) . " of\n$seq\n" if defined $debug;
        return undef;
    }
    return($at + $flanking - 1); # last position of transposon
}

# returns the name of the file
sub BLAT8($$$$$$$$) {
    my ($queriesFile,$dbFile,$tmpdir,$blatcmd,$minScore,$minIdentity,$tileSize,$stepSize) = @_;

    my $blat8File = "$tmpdir/MapTnSeq.$$.".rand() . ".psl";
    # prevent writing to stdout by BLAT
    open(OLD_STDOUT, ">&STDOUT");
    print OLD_STDOUT ""; # prevent warning
    open(STDOUT, ">", "/dev/null");
    my @cmd = ($blatcmd, "-out=blast8", "-t=dna", "-q=dna", "-minScore=$minScore", "-minIdentity=$minIdentity", "-maxIntron=0", "-noTrimA", "-tileSize=$tileSize", "-stepSize=$stepSize", $dbFile, $queriesFile, $blat8File);
    print STDERR join(" ", "Running:", @cmd)."\n" if $debug;
    system(@cmd) == 0 || die "Cannot run $blatcmd: $!";
    close(STDOUT);
    open(STDOUT, ">&OLD_STDOUT") || die "Cannot restore stdout: $!";
    return($blat8File);
}

sub HandleGenomeBLAT($$) {
    my ($rows, $hitsPastEnd) = @_;
    return if scalar(@$rows) == 0;
    my $read = $rows->[0][0];
    die if !defined $read || !exists $nameToBarcode{$read};

    my @hits = ();

    # indexes within besthits entries
    my ($SCAFFOLD,$POSITION,$STRAND,$SCORE,$QBEG,$QEND) = (0,1,2,3,4,5);
    my @besthits = ();

    foreach my $row (@$rows) {
        my ($read2, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $eval, $score) = @$row;
        die unless $read2 eq $read;
        if (scalar(@besthits) == 0 || $score >= $besthits[0][$SCORE] - $delta) {
            # convert from 0-based to 1-based position, and note that sBeg always < sEnd so flip if stranded
            push @besthits,  [ $subject, $sBeg, ($sBeg < $sEnd ? "+" : "-"), $score,
                               $identity, $qBeg, $qEnd ];
            print STDERR "Hit for $read:$qBeg:$qEnd to $subject:$sBeg:$sEnd score $score identity $identity\n"
                if $debug;
        }
    }

    die if @besthits == 0;
    my ($scaffold, $position, $strand, $score, $identity, $qBeg, $qEnd) = @{ $besthits[0] };

    # and output a mapping row (or none)
    if (exists $hitsPastEnd->{$read}) {
        if ($hitsPastEnd->{$read} >= $score - $delta) {
            $nPastEndTrumps++;
            print STDERR "Past end trumps for $read\n" if $debug;
            return;
        } else {
            $nPastEndIgnored++; # ignore weak hit to past-end sequence
            print STDERR "Ignoring hit past end of transposon for $read\n" if $debug;
            $hitsPastEnd->{$read} = 0; # so we do not print out a line for it later on
        }
    }
    # else if no trumping
    print join("\t", $read, $nameToBarcode{$read}, $scaffold, $position, $strand, @besthits == 1 ? 1 : 0,
               $qBeg, $qEnd, $score, $identity)."\n";
    $nMapUnique++ if @besthits == 1;
    $nMapped++;
    return;
}
