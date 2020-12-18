#!/usr/bin/perl -w

#######################################################
## orgGenes.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# This page generates a table of genes for an organism
#
# Required CGI parameters:
# orgId -- which organism to search for

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi = CGI->new;

my $orgId = $cgi->param('orgId');
die "No orgId" if !defined $orgId || $orgId eq "";

my $dbh = Utils::get_dbh();

# gather all of the data needed
my $genes = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE orgId = ? ORDER BY scaffoldId,begin",
                                    {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$genes) == 0;

print "Content-Type:application/x-download\n";
print "Content-Disposition: attachment; filename=organism_${orgId}_genes.tab\n\n";

my @fields = qw{locusId sysName scaffoldId begin end strand type desc};
print join("\t", @fields)."\n";
foreach my $gene (@$genes) {
        my @values = map $gene->{$_}, @fields;
        print join("\t", @values)."\n";
}
