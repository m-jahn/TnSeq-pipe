#!/usr/bin/perl -w
#########################################################################
## StrainFitness.pm -- use db.StrainFitness.orgname and the
## StrainDataSeek table to fetch fitness data for strains in a range of
## locations in the genome.
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Morgan Price
#########################################################################

package StrainFitness;
require Exporter;

use strict;
use DBI;

# Given the path to the cgi_data directory and an open database handle,
# fetch all of the per-strain data for the indicated region of the genome.
# Returns a reference to list of hashes.
sub GetStrainFitness($$$$$$) {
    my ($path, $dbh, $orgId, $scaffoldId, $begin, $end) = @_;
    my $beginKb = int($begin / 1000);
    my $endKb = int(($end + 999) / 1000);

    my ($seek) = $dbh->selectrow_array(qq{ SELECT min(seek) FROM StrainDataSeek
                                           WHERE orgId = ? AND scaffoldId = ? AND kb >= ? AND kb <= ? },
                                       {}, $orgId, $scaffoldId, $beginKb, $endKb);
    return [] if !defined $seek || $seek eq "";
    my $file = "$path/db.StrainFitness.$orgId";
    open(DATA, "<", "$path/db.StrainFitness.$orgId") || die "Cannot read $file";
    my $header = <DATA>;
    chomp $header;
    my @colNames = split /\t/, $header, -1;
    # my %colNames = map {$colNames[$_] => $_} (0..(scalar(@colNames)-1));
    seek(DATA, $seek, 0); # move to absolute position
    my @out = ();
    while(<DATA>) {
        chomp;
        my @F = split /\t/, $_, -1;
        die "Wrong number of columns in\n$_\nfrom $file" unless scalar(@F) == scalar(@colNames);
        my %row = map { $colNames[$_] => $F[$_] } (0..(scalar(@F)-1));
        if ($row{scaffold} ne $scaffoldId || $row{pos} > $end) {
            last;
        } elsif ($row{pos} >= $begin) {
            push @out, \%row;
        }
    }
    close(DATA) || die "Error reading $file";
    return(\@out);
}

1;
