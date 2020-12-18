#!/usr/bin/perl -w
# RunBarSeqLocal.pl -- process a lane of barseq into a table of counts for each strain in each sample
# issuing jobs locally using submitter.pl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../lib";

my $linesPerPiece = 20*1000*1000;
my $minQuality = 0;
my $offset = 0;

my $usage = <<END
Usage:
RunBarSeqLocal.pl [ -indexes BarSeqPrimersH48 | -n25 | -bs3 ]
    organism_directory setname fastq.gz_file_or_directory_with_fastq.gz_files
or
RunBarSeqLocal.pl [ -indexes BarSeqPrimersH48 | -n25 | -bs3 ]
    -in fastq.gz_file_or_directory_with_fastq.gz_files
    -sets organism1_libraryname1_setname1,...,organismN_libraryN_setnameN

Examples:
feba/bin/RunBarSeqLocal.pl g/Keio Keio_ML9_set1 fastq/Keio_ML9_set1
feba/bin/RunBarSeqLocal.pl -indexes feba/primers/BarSeqPrimersH48 g/MR1 MR1_ML3_set1 fastq/MR1_ML3_set1
feba/bin/RunBarSeqLocal.pl -indexes feba/primers/BarSeqPrimersH48 g/psRCH2 psRCH2_ML7_set1 fastq/7341.4.68147.fastq.gz 
feba/bin/RunBarSeqLocal.pl -n25 -in HiSeq_barcodes/FEBA_BS_186 -sets MR1_ML3_set12,Koxy_ML2_set9

RunBarSeq.pl has two phases -- the first phase counts the barcodes in
the fastq file(s).  If the BarSeq was run with the newer style of
primers that have names like IT001, then do not use the -indexes
argument.  In this case it will only accept a directory as input, and
it will assume that each fastq or fastq.gz file in that directory is
already demultiplexed, with a code in the filename such as _IT094_ or
_10_TAGCTT_ to indicate which sample it is from. There can be more
than one file per sample. Usually the files are in the directory, but
if not then it will look across all subdirectories. If using primers
with 2:5 Ns, use the -n25 argument.

(For older experiments, where the primers have names like H01 or M01,
use the -indexes argument to describe which primers were used.
RunBarSeq.pl will use this information to demultiplex the reads.  If
-indexes is used, then RunBarSeq.pl can accept as input either a
directory that contains fastq or fastq.gz files or a single
large_fastq.gz file. Given a single fastq.gz file, it will split it
and process each piece in parallel.)

The second phase aggregates the counts of these barcodes. To save
memory, it ignores barcodes that do not match the pool file, which is
expected to be in

    organism_directory/pool or pool.n10

(If using the -sets argument, organism directory is set to g/organism1, etc.,
and overwriting output files is not allowed. Also you can separate the sets
by ; and optional whitespace instead of by ,.)

Output files are in organism_directory/ --
setname.colsum -- total parsed reads per index
setname.poolcount -- counts for strains in the pool
    based on joining to the pool or pool.n10 file.
setname.codes.ignored -- all the counts for the barcodes
	that did not match the pool.

Other options:
  -minQuality -- Set the minimum quality for each nucleotide in the barcode
	(default is $minQuality)
  -pieceLines $linesPerPiece -- number of lines per piece
  -nosplit -- do not run split, use existing pieces
  -debug -- do not do any work, just show what commands would be run
  -limit -- limit the #reads analyzed per sample (mostly for debugging)
  -offset -- If IT001 is named S97, use -offset 96
  -preseq, -postseq, -nPreExpected -- see MultiCodes.pl
END
    ;

sub maybeRun($); # run command unless $debug is defined
 
my $debug = undef;
{
    my $nosplit = undef;
    my $limitReads = undef;
    my $barcodes = undef;
    my ($n25, $bs3);
    my $setspec = undef;
    my $fastq = undef;
    my ($preseq, $postseq, $nPreExpected);

    die $usage unless GetOptions('debug' => \$debug,
                                 'n25' => \$n25,
                                 'bs3' => \$bs3,
                                 'nosplit' => \$nosplit,
                                 'minQuality=i' => \$minQuality,
                                 'pieceLines=i' => \$linesPerPiece,
				 'limit=i' => \$limitReads,
				 'indexes=s' => \$barcodes,
                                 'preseq=s' => \$preseq,
                                 'postseq=s' => \$postseq,
                                 'nPreExpected=s' => \$nPreExpected,
                                 'in=s' => \$fastq,
                                 'sets=s' => \$setspec,
                                 'offset=i' => \$offset);
    my @gdirs = ();
    my @setnames = ();
    if (defined $setspec) {
      die $usage unless defined $fastq && @ARGV == 0;
      @setnames = split /[,;]\s*/, $setspec;
      my @misc = grep m/^misc_set/i, @setnames;
      print STDERR "Not running combineBarSeq for miscellaneous sets: @misc\n" if @misc > 0;
      @setnames = grep !m/^misc_set/i, @setnames;
      foreach my $setname (@setnames) {
        my @setparts = split /_/, $setname;
        die "Invalid set name: $setname" unless @setparts > 1;
        my $gdir = undef;
        if ($setname =~ m/^(.*)_ML/i) {
          my $org = $1;
          $gdir = "g/$org" if -d "g/$org";
        }
        my $beg = join("_", $setparts[0], $setparts[1]);
        if (!defined $gdir) {
          $gdir = "g/$beg" if -d "g/$beg";
        }
        if (!defined $gdir) {
          # Try to use glob to find the mapping from the first two parts of the library name to an organism
          my @hits = glob("g/*/$beg*");
          die "Could not convert $setname to an organism, and gdir does not exist"
            unless @hits > 0;
          my $hit1 = $hits[0];
          my @parts = split "/", $hit1;
          $gdir = "g/$parts[1]";
          print STDERR "Set $setname matches directory $gdir\n";
        }
        die unless -d $gdir;
        push @gdirs, $gdir;
      }
      # and check that all of these are new
      foreach my $i (0..(scalar(@setnames)-1)) {
        my $setname = $setnames[$i];
        my $gdir = $gdirs[$i];
        my $out = "$gdir/$setname";
        foreach my $file ("$out.colsum", "$out.poolcount") {
          die join("\n",
                   "Error: output for $gdir $setname already exists!",
                   "See file $file",
                   "Please remove the colsum and poolcount files or correct the set name",
                   "")
            if -s $file; # if not empty
        }
      }
    } else {
      die $usage unless @ARGV == 3;
      my ($gdir, $setname);
      ($gdir, $setname, $fastq) = @ARGV;
      die "No such directory: $gdir" unless -d $gdir;
      push @setnames, $setname;
      push @gdirs, $gdir;
    }
    die "No such file: $fastq" unless -e $fastq;
    die "No such file: $barcodes" unless !defined $barcodes || -e $barcodes;
    # And check that all gdir are unique
    my %gdirs = map { $_ => 1 } @gdirs;
    die "Each organism can be included only once in the -sets argument\n"
      unless scalar(keys %gdirs) == scalar(@gdirs);

    my @poolfiles = ();
    foreach my $gdir (@gdirs) {
      my $poolfile = "$gdir/pool";
      $poolfile = "$poolfile.n10" if !-e $poolfile;
      die "Cannot find $poolfile (or withut .n10)" unless -e $poolfile;
      push @poolfiles, $poolfile;
    }

    # design the parts
    my $prefix = $fastq;
    my @parts;
    my @codes;
    my $codeGlob;
    my $setname1 = $setnames[0]; # for naming the parts
    if (-d $fastq) {
	@parts = glob("$fastq/*fastq.gz");
	@parts = grep { !m/^[.]/ }  @parts;
        # Ignored Undetermined_* files before deciding if subdirectories need to be searched
        @parts = grep { !m/Undetermined_/ } @parts;
	if (scalar(@parts) == 0) {
	    die "No *.gz files in $fastq" if defined $barcodes;
	    # Indexed runs sometimes have one directory per sample, with fastq.gz file(s) within each directory
	    @parts = glob("$fastq/*/*.fastq.gz");
	    @parts = grep !m"[.]/", @parts; # remove hidden files
	    die "Cannot find *.gz files in $fastq or its subdirectories" unless scalar(@parts) > 0;
	    $codeGlob = "$fastq/*/*.codes";
	} else {
	    $codeGlob = "$fastq/*.codes";
	}
	@codes = map { my $i = $_; $i =~ s/.fastq.gz$/.codes/; $i; } @parts;
	@codes = grep { -e $_ } @codes; # only the files that were previously made
    } else {
	die "Cannot demultiplex the file $fastq without the -indexes option" unless defined $barcodes;
	$prefix = $fastq; $prefix =~ s/fastq[.]gz$//;
	$prefix = "$prefix.parts";
	if (! -d $prefix) {
	    mkdir($prefix) || die "Cannot mkdir $prefix; $!";
	}
        @parts = glob("$prefix/${setname1}_BarSeq.part*[0-9]");
        $codeGlob = "$prefix/${setname1}_BarSeq.part*[0-9].codes";
        @codes = glob($codeGlob); # only the preexisting files
    }
    print STDERR "See " . scalar(@parts) . " parts and " . scalar(@codes) . " codes files; codes files will be overwritten\n"
	if scalar(@codes) > 0;

    # splitting if necessary
    if (! -d $fastq) {
        if (defined $nosplit && scalar(@parts) > 0) {
            print STDERR "Skipping the split step, just redoing the analysis with " . scalar(@parts) . " existing pieces\n";
            maybeRun("rm $codeGlob")
                if scalar(@codes) > 0 && ! defined $debug;
        } else {
            # do the splitting
            maybeRun("rm $prefix/${setname1}_BarSeq.part*") if (scalar(@parts) > 0 || scalar(@codes) > 0);
            my $cmd = "gunzip -c $fastq | split -l $linesPerPiece -d - $prefix/${setname1}_BarSeq.part";
            maybeRun($cmd);
            system("touch $prefix/${setname1}_BarSeq.part00") if defined $debug;
	    @parts = glob("$prefix/${setname1}_BarSeq.part*[0-9]"); # now @parts is actually there
        }
    }

    # build the list of commands, and update @codes to be what we'll make
    @codes = (); # will update with expected results
    my $cmdsfile = "$prefix/${setname1}_BarSeq.codecmds";
    maybeRun("rm $cmdsfile*") if -e $cmdsfile;
    maybeRun("rm -f $codeGlob");
    open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
    foreach my $i (@parts) {
	print STDERR "Considering part $i\n" if $debug;

        my $corecmd = "$Bin/MultiCodes.pl -minQuality $minQuality";
        $corecmd .= " -n25" if defined $n25;
        $corecmd .= " -bs3" if defined $bs3;
	$corecmd .= " -limit $limitReads" if defined $limitReads;
        $corecmd .= " -preseq $preseq" if defined $preseq;
        $corecmd .= " -postseq $postseq" if defined $postseq;
        $corecmd .= " -nPreExpected $nPreExpected" if defined $nPreExpected;
	if (defined $barcodes) {
	    $corecmd .= " -primers $barcodes";
	} else {
	    my @path = split "/", $i;
	    my $name = pop @path;
	    my @pieces = split /[._]/, $name;
	    my @indexes = grep m/^IT\d+$/, @pieces;
	    my $index = undef;
            my $sampleNum = undef;
	    # sometimes ITO (capital O) not IT0 (numeral zero)
	    my @indexes2 = grep m/^ITO\d+$/, @pieces;
	    my @indexes3 = grep m/IT0\d\d/, @pieces;

	    if (@indexes == 1) {
              $index = $indexes[0];
	    } elsif (@indexes2 == 1){
              # e.g. HiSeq_barcodes/FEBA_BS_81/Sample_FEBA_BS_81_ITO88/FEBA_BS_81_ITO88_GGCCTG_L001_R1_001.fastq.gz
              $index = $indexes2[0];
              $index =~ s/ITO/IT0/;
            } elsif (@indexes3 == 1) {
              $indexes3[0] =~ m/(IT0\d\d)/ || die;
              $index = $1;
	    } elsif ($name =~ m/_(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/
		     || $name =~ m/_Index(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/i
                     || $name =~ m/_Index(\d+)_S\d+_/
                     || $i =~ m!_IT(\d\d\d)[/_.]!) {
              # e.g. FEBA_BS_60_10_TAGCTT_L001_R1_001.fastq.gz
              # e.g. FEBA_BS_125_Index10_S10_L001_R1_001.fastq.gz
              # e.g. FEBA_BS_195_IT001/FEBA_BS_195_S97_L002_R1_001.fastq.gz
              $sampleNum = $1;
            } elsif ($name =~ m/_(\d+)_S\d+_L\d+_/) {
              # e.g. FEBA_BS_117_24_S120_L002_R1_001.fastq.gz is IT024
              $sampleNum = $1;
            } elsif ($name =~ m/^(\d+)_S\d+_L\d+_R\d+_/) {
              # e.g. 18_S1_L001_R1_001.fastq.gz is IT018
              $sampleNum = $1;
            } elsif ($name =~ m/^[A-Z]*\d*_S(\d+)_*L*\d*_R\d+/) {
              # e.g. A10_S106_L003_R1_001.fastq.gz
              # e.g. A10_S10_R1_001.fastq.gz
              $sampleNum = $1;
	    } elsif ($name =~ m/undetermined/i) {
              print STDERR "Skipping $name\n";
              next;
	    } else {
              die "Cannot identify the index ITnnn from file $i";
	    }
            if (defined $sampleNum && !defined $index) {
              die "Invalid sample number $sampleNum with offset $offset from file $name\n"
                unless $sampleNum - $offset >= 1;
              $index = sprintf("IT%03d", $sampleNum - $offset);
            }
            die unless defined $index;
	    $corecmd .= " -index $index";
	}
        if ($i =~ m/[.]gz$/) {
            my $out = $i;
            $out =~ s/.fastq.gz$//;
	    push @codes, "$out.codes";
            print CMDS "zcat $i | $corecmd -out $out >& $out.log"."\n";
        } else {
	    push @codes, "$i.codes";
            print CMDS "$corecmd -out $i < $i >& $i.log"."\n";
        }
    }
    close(CMDS) || die "Error writing to $cmdsfile";

    # run them all
    maybeRun("$Bin/submitter.pl $cmdsfile");

    # combine the results
    if (!defined $debug) {
	# Parse all the log files
	my @logs = map { my $out = $_; $out =~ s/[.]codes$/.log/; $out; } @codes;
	my $nReads = 0;
	my $nMulti = 0;
	my $nFiles = 0;
	my $nUsable = 0;
	foreach my $log (@logs) {
	    if (-e $log) {
		$nFiles++;
		open(LOG, "<", $log) || die "Cannot read $log";
		my $nLines = 0; # numbers of lines parsed with counts of reads
		while (<LOG>) {
		    if (m/^Reads\s+(\d+)\s+Multiplexed\s+(\d+)\s+Usable\S*\s+(\d+)\b/) {
			$nReads += $1;
			$nMulti += $2;
			$nUsable += $3;
			$nLines++;
		    }
		}
		close(LOG) || die "Error reading $log";
		print STDERR "Warning: no Reads/Multiplexed/Usable line in $log\n"
		    if $nLines == 0;
	    } else {
		print STDERR "Warning: no log file $log\n";
	    }
	}
	print STDERR sprintf("Total reads (in millions) %.1f Multi %.1f (%.1f%%) Usable %.1f (%.1f%%) from %d files%s\n",
			     $nReads/1e6,
			     $nMulti/1e6, 100.0*$nMulti/($nReads+0.1),
			     $nUsable/1e6, 100.0*$nUsable/($nReads+0.1),
			     $nFiles,
			     defined $setspec ? "" : "for $setname1");
    }
    foreach my $i (0..(scalar(@gdirs)-1)) {
      my $path = $gdirs[$i] . "/" . $setnames[$i];
      print STDERR "Combining codes files $codeGlob into $path\n";
      my $cmd = "$Bin/combineBarSeq.pl $path $poolfiles[$i] $codeGlob";
      maybeRun($cmd);
    }
}

sub maybeRun($) {
    my ($cmd) = @_;
    if (defined $debug) {
        print STDERR "Would run: $cmd\n";
    } else {
        system($cmd) == 0 || die "script failed: $cmd";
    }
}
