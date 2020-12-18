#!/usr/bin/perl -w

#######################################################
## createExpData.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page generates output to STDOUT and presents it as
# a download, which includes all of the experimental data
# for a single organism. It is run from org.cgi. 
# 
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
$orgId = "" if !defined $orgId;

my $dbh = Utils::get_dbh();

# gather all of the data needed

my $gene = $dbh->selectall_arrayref("SELECT orgId, locusId, sysName, gene FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

my $exp = $dbh->selectall_arrayref("SELECT * FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);



print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=exp_organism_$orgId.txt\n\n";

# open my $logFile, '>', "organism_$orgId.txt" or die "error trying to (over)write: $!";

# print the header row
my @fields = qw{orgId expName expDesc timeZeroSet
                num nMapped nPastEnd nGenic nUsed gMed gMedt0 gMean cor12 mad12 mad12c mad12c_t0 opcor adjcor gccor maxFit
                expGroup expDescLong mutantLibrary person dateStarted setName seqindex
                media temperature pH vessel aerobic liquid shaking
                condition_1 units_1 concentration_1
                condition_2 units_2 concentration_2
                condition_3 units_3 concentration_3
                condition_4 units_4 concentration_4
                growthPlate growthWells};
print join("\t", @fields)."\n";

# print the data row by row
foreach my $row(@$exp) {
  my @out = ();
  foreach my $field (@fields) {
    die "Invalid field $field" unless exists $row->{$field};
    push @out, $row->{$field};
  }
  print join("\t", @out)."\n";
}
