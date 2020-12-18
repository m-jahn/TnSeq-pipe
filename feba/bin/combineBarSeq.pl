#!/usr/bin/perl -w
use strict;
# Given a bunch of files with counts for barcodes, pick out the
# ones from the pool and make a small table of their counts (out.poolcount).
# Also makes a table of total counts (out.colsum).
# Expects either that all columns are in all files,
# or that each file has just one column

use POSIX; # for floor, ceil
sub median;
sub ShortList; # prettier form of a list of indexes

my $save_ignore = 0; # make a file of ignored lines (out.codes.ignored) ?

{
    die "Usage: combineBarSeq.pl out pool_file codesfiles...\n"
      . "  Codes files may be gzipped.\n"
      . "  Writes out.poolcount and out.colsum.\n"
	unless @ARGV >= 3;
    my $out = shift @ARGV;
    my $poolfile = shift @ARGV;
    my @codesFiles = @ARGV;

    my %pool = (); #  rcbarcode to barcode,scaffold,strand,pos
    open(POOL, "<", $poolfile) || die "Cannot read $poolfile";
    while(<POOL>) {
	chomp;
	my ($barcode,$rcbarcode,undef,undef,$scaffold,$strand,$pos) = split /\t/, $_;
	next if $barcode eq "barcode";
	die "Invalid barcode $barcode" unless $barcode =~ m/^[ACGT]+$/;
	die "Invalid rcbarcode $rcbarcode" unless $rcbarcode =~ m/^[ACGT]+$/;
	die "Invalid position $pos" unless $pos eq "" || $pos =~ m/^\d+$/;
	die "Invalid strain $strand" unless $strand eq "+" || $strand eq "-" || $strand eq "";
	die "Duplicate rcbarcode $rcbarcode" if exists $pool{$rcbarcode};
	$pool{$rcbarcode} = [$barcode,$scaffold,$strand,$pos];
    }
    close(POOL) || die "Error reading $poolfile";
    die "No entries in pool file $poolfile" unless scalar(keys %pool) > 0;

    if ($save_ignore) {
	open(IGNORE, ">", "$out.codes.ignored") || die "Cannot write to $out.codes.ignored";
    }

    my %counts = (); # rcbarcode to vector of counts
    my @indexes = (); # vector of names of samples
    my @colSums = (); # samples to total number of counts
    my @colSumsUsed = (); # samples to total number of counts for used barcodes
    my $nSamples = 0;
    my $oneperfile = 0;
    my $thisIndex = undef; # if in one-per-file mode, which index is this file?

    my $nUsed = 0;
    my $nIgnore = 0;
    foreach my $file (@codesFiles) {
      if ($file =~ m/[.]gz$/) {
        open(IN, "zcat $file |") || die "Cannot zcat $file";
      } else {
	open(IN, "<", $file) || die "Cannot read $file";
      }
      my $header = <IN>;
      chomp $header;
      my @cols = split /\t/, $header;
      my $first = shift @cols;
      die "Not a barcode counts file -- $file" unless $first eq "barcode";
      die "No index columns in $file" unless @cols > 0;
      if ($nSamples == 0) { # first file
        @indexes = @cols;
        $nSamples = scalar(@indexes);
        @colSums = (0) x $nSamples;
        @colSumsUsed = (0) x $nSamples;
        if (scalar(@cols) == 1) {
          $oneperfile = 1;
          $thisIndex = 0;
        }
      } elsif ($oneperfile) {
        die "More than one data column in $file" unless scalar(@cols) == 1;
        my %oldcol = map { $indexes[$_] => $_ } (0..(scalar(@indexes)-1));
        my $thisCol = $cols[0];
        if (exists $oldcol{ $thisCol }) {
          $thisIndex = $oldcol{ $thisCol };
        } else {
          $nSamples++;
          push @indexes, $thisCol;
          $thisIndex = scalar(@indexes)-1;
          push @colSums, 0;
          push @colSumsUsed, 0;
        }
      } else {
        die "Wrong number of columns in $file" unless scalar(@cols) == scalar(@indexes);
        foreach my $i (0..(scalar(@cols)-1)) {
          die "Index mismatch in $file vs. $codesFiles[0] -- $cols[$i] vs. $indexes[$i]"
            unless $cols[$i] eq $indexes[$i];
        }
      }
      my $nThisFile = 0;

      while(<IN>) {
        $nThisFile++;
        chomp;
        my @F = split /\t/, $_;
        my $barcode = shift @F; # actually rcbarcode
        # note am allowing N in barcode but not in pool
        die "Invalid barcode: $barcode" unless $barcode =~ m/^[ACGTN]+$/; 
        die "Wrong number of columns in $file" unless scalar(@F) == ($oneperfile ? 1 : $nSamples);
        if (exists $pool{$barcode}) {
          $nUsed++;
          if ($oneperfile) {
            $colSumsUsed[$thisIndex] += $F[0];
            $counts{$barcode}[$thisIndex] += $F[0];
          } else {
            for (my $i = 0; $i < $nSamples; $i++) {
              $colSumsUsed[$i] += $F[$i];
            }
            if (exists $counts{$barcode}) {
              my $row = $counts{$barcode};
              for (my $i = 0; $i < $nSamples; $i++) {
                $row->[$i] += $F[$i];
              }
            } else {
              $counts{$barcode} = \@F;
            }
          }
        } else {
          print IGNORE join("\t",$barcode,@F)."\n" if $save_ignore;
          $nIgnore++;
        }
        if ($oneperfile) {
          $colSums[$thisIndex] += $F[0];
        } else {
          for (my $i = 0; $i < $nSamples; $i++) {
            $colSums[$i] += $F[$i];
          }
        }
      }
      close(IN) || die "Error reading from $file";
      print STDERR "Warning: no entries in $file\n" if $nThisFile == 0;
    }
    print STDERR sprintf("Pool %s entries %.1fK saw %.1fK lines ignored %.1fK lines from %d files\n",
                         $poolfile, scalar(keys %pool)/1000.0, $nUsed/1000.0, $nIgnore/1000.0, scalar(@codesFiles));
    if ($save_ignore) {
      close(IGNORE) || die "Error writing to $out.codes.ignored";
    }

    # categorize reads
    my @lowcountI = (); # indexes with < 200,000 reads
    my @lowcounts = (); # those actual counts
    my @fractions = (); # fraction used for indexes with plenty of reads only
    my @fractionsS = (); # fraction used for successful indexes (enough reads and f > 0.25)
    my @lowhitI = (); # indexes with fraction < 0.25
    my @okI = (); # indexes with enough reads and fraction above 0.25
    my $totalReads = 0;
    foreach my $i (0..(scalar(@indexes)-1)) {
	$totalReads += $colSums[$i];
	if ($colSums[$i] < 200*1000) {
	    push @lowcountI, $indexes[$i];
	    push @lowcounts, $colSums[$i];
	} else {
	    my $fraction = $colSumsUsed[$i] / $colSums[$i];
	    push @fractions, $fraction;
	    if ($fraction < 0.25) {
		push @lowhitI, $indexes[$i];
	    } else {
		push @okI, $indexes[$i];
                push @fractionsS, $fraction;
	    }
	}
    }
    print STDERR sprintf("  Indexes %d Success %d LowCount %d LowFraction %d Total Reads (millions): %.3f \n",
			 scalar(@indexes), scalar(@okI), scalar(@lowcountI), scalar(@lowhitI), $totalReads/1e6);
    print STDERR sprintf("  Median Fraction (LowCount Excluded) %.3f Median for Success %.3f\n",
			 &median(@fractions),
                         scalar(@fractionsS) > 0 ? &median(@fractionsS) : "NaN")
      if scalar(@fractions) > 0;
    print STDERR "  Success " . &ShortList(@okI) . "\n" if @okI > 0;
    print STDERR sprintf("  LowCount (median %d) %s\n",
			 &median(@lowcounts), &ShortList(@lowcountI))
	if @lowcountI > 0;

    open(COUNT, ">", "$out.poolcount") || die "Cannot write to $out.poolcount";
    print COUNT join("\t", "barcode", "rcbarcode", "scaffold", "strand", "pos", @indexes)."\n";
    # Note -- BarSeqR.pl expects the ordering to be the same across runs.
    # Older versions of the code used the perl hash order, which was risky.
    # In perl 5.18 hash ordering is randomized and this broke.
    # So, now we explicitly sort the barcodes
    # However, if you analyze results from old and new combineBarSeq.pl runs
    # then they will collide.
    my @sorted = sort { my ($bc1,$sc1,$strand1,$pos1) = @{ $pool{$a} };
                        my ($bc2,$sc2,$strand2,$pos2) = @{ $pool{$b} };
                        $sc1 cmp $sc2
                            || ($sc1 eq "pastEnd" ? 0 : $pos1 <=> $pos2)
                            || $strand1 cmp $strand2
                            || $bc1 cmp $bc2;
    } keys %pool;
    foreach my $rcbarcode (@sorted) {
	my ($barcode,$scaffold,$strand,$pos) = @{ $pool{$rcbarcode} };
	my $counts = $counts{$rcbarcode};
	my @out = ($barcode,$rcbarcode,$scaffold,$strand,$pos);
	if (defined $counts) {
	    for (my $i = 0; $i < $nSamples; $i++) {
		$counts->[$i] += 0; # ensure vector is full length, not undef;
	    }
	    push @out, @$counts;
	} else {
	    push @out, (0) x $nSamples;
	}
	print COUNT join("\t", @out)."\n";
    }
    close(COUNT) || die "Error writing to $out.poolcount";
    print STDERR "  Wrote $out.poolcount\n";

    open(SUM, ">", "$out.colsum") || die "Cannot write to $out.colsum";
    print SUM join("\t", qw{Index nReads nUsed fraction})."\n";
    foreach my $i (0..(scalar(@indexes)-1)) {
	print SUM join("\t", $indexes[$i], $colSums[$i], $colSumsUsed[$i], $colSumsUsed[$i]/(0.5 + $colSums[$i]))."\n";
    }
    close(SUM) || die "Error writing to $out.colsum";
    print STDERR "  Wrote $out.colsum\n";
}

sub median {
    return undef if scalar(@_) == 0;
    my @sorted = sort { $a <=> $b } @_;
    my $midindex = (scalar(@sorted) - 1)/2;
    return ($sorted[floor($midindex)] + $sorted[ceil($midindex)]) / 2.0;
}

# assume indexes are of the form alphabetic_prefix followed by a number
# shortens it by replacing prefix1 prefix2 prefix3 with prefix1:3
# outputs a pretty string to print
sub ShortList {
    my @list = @_;
    @list = sort(@list);
    return join(" ", @list) unless $list[0] =~ m/^([a-zA-Z]+)\d+$/;
    my $prefix = $1;
    my $prelen = length($prefix);
    my @numbers = map { substr($_, 0, $prelen) eq $prefix && substr($_, $prelen) =~ m/^\d+$/ ?
			    substr($_, $prelen) : undef } @list;
    my $lastno = undef;
    my $inrun = 0;
    my $sofar = "";
    foreach my $i (0..(scalar(@list)-1)) {
	my $string = $list[$i];
	my $number = $numbers[$i];
	if (defined($number) && defined($lastno) && $number == $lastno+1) {
	    $sofar .= ":$number" if $i == scalar(@list)-1; # terminate if at end of list
	    $inrun = 1; # start or continue a run
	    $lastno = $number;
	} else {
	    if (defined $lastno) {
		$sofar .= ":$lastno" if $inrun;
		$lastno = undef;
	    }
	    $inrun = 0;
	    $sofar .= " " unless $sofar eq "";
	    $sofar .= $string;
	    $lastno = $number;
	}
    }
    return($sofar);
}
