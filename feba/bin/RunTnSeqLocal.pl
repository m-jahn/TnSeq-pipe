#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Which;

my $usage = <<END
Usage: RunTnSeqLocal.pl [ -nosplit ] [ -limit nMaxReadsPerFile ] [ -blat blatcmd ]
    nickname library_name modelsuffix fastq_directory_or_gz
Must be run in the parent directory of g/nickname
The file model_modelsuffix must exist in the feba primers/ directory.
Writes to g/nickname/library_name
-limit and -blat are passed on to MapTnSeq.pl
END
    ;

{
    my $limit = undef;
    my $blat = -e "$Bin/blat" ? "$Bin/blat" : "blat";
    my $nosplit = undef;
    die $usage unless GetOptions('limit=i' => \$limit,
				 'blat=s' => \$blat,
	                         'nosplit' => \$nosplit)
	&& @ARGV==4;
    my ($nickname, $library, $modelsuffix, $fastqdir) = @ARGV;

    die "blat is not on the path!" if ($blat eq "blat" && !defined which("blat"));

    die "No such directory: g/$nickname" unless -d "g/$nickname";
    my $modelfile = "$Bin/../primers/model_".$modelsuffix;
    die "No such file: $modelfile" unless -e $modelfile;

    die "Not such file: $fastqdir" unless -e $fastqdir;
    if (! -d $fastqdir) {
	# split the gz file into pieces
	my $file = $fastqdir;
	my $gzipped = 0;
	if ($file =~ m/[.]gz$/) {
	    $gzipped = 1;
	}
	$fastqdir =~ s/fastq[.]gz$//;
	$fastqdir =~ s/fastq$//;
	$fastqdir .= ".parts";
	if (! -d $fastqdir) {
	    mkdir($fastqdir) || die "Cannot mkdir $fastqdir: $!";
	    undef $nosplit;
	}
	my $zcat_cmd = $gzipped ? "zcat $file" : "cat $file";
	my $num_lines = 20*1000*1000;
	if (defined $nosplit ) {
	    print STDERR "Skipping splitting $file into pieces in $fastqdir\n";
	} else {
	    print STDERR "Splitting $file into pieces in $fastqdir\n";
	    system("$zcat_cmd | split --lines=$num_lines --numeric-suffixes - $fastqdir/TnSeq_part.fastq") == 0
		|| die "split failed";
	}
    }
    my $cmdsfile = "$fastqdir/RunTnSeqLocal.cmds";
    system("rm ${cmdsfile}*") if -e $cmdsfile;

    my @inputs = ();
    @inputs = glob("$fastqdir/*.fastq.gz");
    @inputs = glob("$fastqdir/*.fastq") if (scalar(@inputs) == 0);
    @inputs = glob("$fastqdir/TnSeq_part.fastq*[0-9]") if (scalar(@inputs) == 0);
    die "Cannot find input fastq or fastq.gz files in $fastqdir" if scalar(@inputs) == 0;
    my @mapped = ();

    open(CMDS, ">", $cmdsfile) || die "Error writing to $cmdsfile";
    foreach my $file (@inputs) {
	my $base = $file;
	if ($file =~ m/[.]fastq[.]gz/) {
	    $base =~ s/[.]fastq[.]gz$//;
	} elsif ($file =~ m/[.]fastq[.]/) {
	    $base =~ s/[.]fastq//;
	} elsif ($file =~ m/part[.]fastq\d+/) {
	    $base =~ s/part[.]fastq/part/;
	}
	push @mapped, "$base.mapped";
	my $cmd = "$Bin/MapTnSeq.pl -genome g/$nickname/genome.fna -model $modelfile -first $file";
	$cmd .= " -limit $limit" if defined $limit;
	$cmd .= " -blat $blat" if defined $blat;
	$cmd .= " > $base.mapped.tmp && mv $base.mapped.tmp $base.mapped";
	unlink("$base.mapped") if -e "$base.mapped";
	print CMDS "$cmd\n";
    }
    close(CMDS) || die "Error writing to $cmdsfile";
    print STDERR "Wrote " . scalar(@inputs) . " mapping jobs to $cmdsfile\n";
    system("$Bin/submitter.pl",$cmdsfile) == 0
	|| die "Error running $Bin/submitter.pl on $cmdsfile\n";

    system("(grep -h Reads $cmdsfile-*.log; grep -h Proportion $cmdsfile-*.log) >& g/$nickname/$library.reads");
    print STDERR "Wrote g/$nickname/$library.reads\n";

    my $cmd = "$Bin/DesignRandomPool.pl -minN 10 -genes g/$nickname/genes.tab -pool g/$nickname/$library.pool "
	. join(" ",@mapped) . " >& g/$nickname/$library.pool.stats";
    print STDERR "Running: $cmd\n";
    system($cmd) == 0
	|| die "DesignRandomPool.pl failed";
    print STDERR "Wrote g/$nickname/$library.pool and .pool.stats\n";
}
