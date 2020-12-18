#!/usr/bin/perl -w
use strict;

my $usage = <<END
Usage: resortPoolCount.pl filename1 ... filenameN

This script updates setname.poolcount files from before Feb 2016 to be
compatible with the new biological order. Saves the original file in
setname.poolcount.oldsort and issues a warning and does NOT overwrite if
that file already exists.
END
    ;

die $usage if @ARGV == 0 || $ARGV[0] =~ m/^-/;
my @files = @ARGV;

foreach my $file (@files) {
    die "No such file: $file\n" unless -e $file;
}
foreach my $file (@files) {
    my $tmp = "$file.tmp.$$";
    my $old = "$file.oldsort";
    if (-e $old) {
        print STDERR "$file.old already exists, has this already been resorted?\n Not updating $file\n";
        next;
    }
    #else
    open(FILE, "<", $file) || die "Cannot read $file";
    my $header = <FILE>;
    chomp $header;
    my @colnames = split /\t/, $header, -1;
    my @expected = qw{barcode rcbarcode scaffold strand pos};
    foreach my $i (0..(scalar(@expected)-1)) {
        die "Do not see expected header name $expected[$i] at position $i in $file"
            unless $colnames[$i] eq $expected[$i];
    }
    my @rows = ();
    while(<FILE>) {
        chomp;
        my @F = split /\t/, $_, -1;
        die "Wrong number of columns in $file" unless scalar(@F) == scalar(@colnames);
        push @rows, \@F;
    }
    close(FILE) || die "Error reading $file";
    @rows = sort {
        my ($bc1,undef,$sc1,$strand1,$pos1) = @$a;
        my ($bc2,undef,$sc2,$strand2,$pos2) = @$b;
        $sc1 cmp $sc2
            || ($sc1 eq "pastEnd" ? 0 : $pos1 <=> $pos2)
            || $strand1 cmp $strand2
            || $bc1 cmp $bc2;
    } @rows;
    open(TMP, ">", $tmp) || die "Cannot write to $tmp";
    print TMP join("\t", @colnames)."\n";
    foreach my $row (@rows) {
        print TMP join("\t", @$row)."\n";
    }
    close(TMP) || die "Error writing to $tmp";
    rename($file,$old) || die "Cannot rename $file to $old";
    rename($tmp,$file) || die "Cannot rename $tmp to $file";
    print STDERR "Wrote $file and saved old sort in $old\n";
}
print STDERR "Done\n";

