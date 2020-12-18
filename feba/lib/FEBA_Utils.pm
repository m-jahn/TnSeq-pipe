
# FEBA_Utils.pm -- Utilities for the command-line scripts.
# Copyright (c) 2017 All Rights Reserved by UC Berkeley
# Authors: Morgan Price at LBL

package FEBA_Utils;
require Exporter;
use strict;
use File::stat;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( ReadTable ReadColumnNames ReadFasta NewerThan ReadFastaDesc ReadFastaEntry reverseComplement
              sqlite_quote );

# filename and list of required fields => list of hashes, each with field->value
# The list can be a single name or a reference to a list
# (It used to be an actual list; not sure how that went wrong with perl prototypes?)
sub ReadTable($*) {
    my ($filename,@required) = @_;
    if (scalar(@required) == 1 && ref $required[0]) {
        @required = @{ $required[0] };
    }
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $headerLine = <IN>;
    $headerLine =~ s/[\r\n]+$//; # for DOS
    # Check for Mac style files -- these are not supported, but give a useful error
    die "Tab-delimited input file $filename is a Mac-style text file, which is not supported\n"
        . "Use\ndos2unix -c mac $filename\nto convert it to a Unix text file.\n"
        if $headerLine =~ m/\t/ && $headerLine =~ m/\r/;
    my @cols = split /\t/, $headerLine, -1;
    my %cols = map { $cols[$_] => $_ } (0..(scalar(@cols)-1));
    foreach my $field (@required) {
	die "No field $field in $filename" unless exists $cols{$field};
    }
    my @rows = ();
    while(my $line = <IN>) {
	$line =~ s/[\r\n]+$//;
	my @F = split /\t/, $line, -1;
	die "Wrong number of columns in:\n$line\nin $filename"
	    unless scalar(@F) == scalar(@cols);
	my %row = map { $cols[$_] => $F[$_] } (0..(scalar(@cols)-1));
	push @rows, \%row;
    }
    close(IN) || die "Error reading $filename";
    return @rows;
}

# filename to list of column names
sub ReadColumnNames($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $line = <IN>;
    close(IN) || die "Error reading $filename";

    $line =~ s/[\r\n]+$//; # for DOS
    my @cols = split /\t/, $line;
    return @cols;
}

# returns a reference to a hash of name => sequence
sub ReadFasta ($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my %seqs = ();
    my $name = undef;
    while(<IN>) {
	chomp;
	if (m/>(\S+)/) {
	    die "Duplicate entry for $1 in filename" if exists $seqs{$1};
	    $name = $1;
	} else {
	    die "Unexpected line $_" unless defined $name;
	    $seqs{$name} .= $_;
	}
    }
    close(IN) || die "Error reading $filename";
    return(\%seqs);
}

# Returns a hash containing either "error"
# or hashes of "seq", "desc", and "len"
sub ReadFastaDesc($) {
    my ($file) = @_;
    my %seq = ();
    my %desc = ();
    my $name = undef;
    open(FAA, "<", $file) || return('error' => "Cannot read $file" );
    while(<FAA>) {
        s/[\r\n]+$//;
        if (m/^>(.*)$/) {
            my $header = $1;
            if ($header =~ m/^(\S+)\s+(\S.*)$/) {
                $name = $1;
                $desc{$name} = $2;
            } else {
                return('error' => "bad header for sequence:\n$header\n") unless $header =~ m/^\S+$/;
                $name = $header;
                $desc{$name} = $header;
            }
            return('error' => "Duplicate sequence id:\n$name\n") if exists $seq{$name};
            $seq{$name} = "";
        } else {
            return('error' => "sequence before header:\n$_\n") unless defined $name;
            s/\s//g;
            $seq{$name} .= $_;
        }
    }
    close(FAA) || return('error' => "Error reading $file");
    my %len = ();
    while (my ($name,$seq) = each %seq) {
        $len{$name} = length($seq);
        return('error' => "No sequence for id:\n$name\n") if ($len{$name} == 0);
    }
    return("desc" => \%desc, "seq" => \%seq, "len" => \%len);
}

# Read one entry at a time from a fasta file
# The first argument is a hash to keep track of saved state, i.e.:
#   my $state = {};
#   while(my ($header,$sequence) = ReadFastaEntry($fh,$state)) { ... }
# (header will have the ">" removed)
sub ReadFastaEntry {
  my ($fh, $state) = @_;
  die unless ref $state;
  return () if exists $state->{DONE}; # end state
  # initialization
  if (!defined $state->{header}) {
    $state->{header} = "";
    $state->{sequence} = "";
  }
  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ m/^>(.*)/) {
      my $old_header = $state->{"header"};
      my $old_sequence = $state->{"sequence"};
      $state->{"header"} = $1;
      die "Empty header in $line" if $state->{header} eq "";
      $state->{"sequence"} = "";
      return ($old_header, $old_sequence) if $old_header ne "";
    } else {
      die "Unexpected sequence with no header" if $state->{"header"} eq "";
      $line =~ s/ //g;
      $line = uc($line);
      # allow - or . as used in alignments and * as used for stop codons
      die "Invalid sequence line $line" unless $line =~ m/^[A-Z*.-]*$/;
      $state->{sequence} .= $line;
    }
  }
  # reached end of file
  $state->{DONE} = 1;
  return () if $state->{header} eq ""; # empty file
  return ($state->{header}, $state->{sequence});
}

sub NewerThan($$) {
    my ($file1, $file2) = @_;
    die "Invalid arguments to NewerThan" unless defined $file1 && defined $file2;
    die "No such file: $file2" unless -e $file2;
    return 0 unless -e $file1;
    return stat($file1)->mtime >= stat($file2)->mtime ? 1 : 0;
}

sub reverseComplement($) 
{
    my $seq = shift;
    chomp $seq;
        my $origSeq=$seq;

    die "Invalid sequence \"$origSeq\" in reverseComplement" if ( ($seq =~ 
tr/RYKMSWBDHVNATCGXrykmswbdhvnatcg-/YRMKWSVHDBNTAGCXyrmkwsvhdbntagc-/) != 
length($seq) );
    $seq = reverse $seq;
    return $seq;
}

# This is for uploading to sqlite3. 
# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub sqlite_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}

1;
