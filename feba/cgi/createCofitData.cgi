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
# a download, which includes all of the experimental FITNESS data
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

my $genes = $dbh->selectall_hashref("SELECT * FROM Gene WHERE orgId = ?",
                                    "locusId", # field to index hash by
                                    {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(keys %$genes) == 0;
my $cofitness = $dbh->selectall_arrayref(qq{ SELECT * FROM Cofit WHERE orgId = ? AND rank <= 20
                                             ORDER BY locusId,rank },
                                        {Slice => {}}, $orgId);
my $cons = $dbh->selectall_arrayref(qq{ SELECT DISTINCT locusId, hitId FROM ConservedCofit
                                        WHERE orgId = ? },
                                    {}, $orgId);
my %cons = (); # locusId => hitId > 1 if conserved
foreach my $row (@$cons) {
    $cons{ $row->[0] }{ $row->[1] } = 1;
}

print "Content-Type:application/x-download\n";
print "Content-Disposition: attachment; filename=cofit_organism_$orgId.txt\n\n";

# print the header row
print join("\t", qw{orgId locusId sysName name desc hitId hitSysName hitName hitDesc rank cofit conserved})."\n";

# print the data row by row
foreach my $row (@$cofitness) {
    my $locusId = $row->{locusId};
    my $gene = $genes->{$locusId};
    die "Unknown locusId $locusId in Cofitness.locusId" if !defined $gene;
    my $hitId = $row->{hitId};
    my $hit = $genes->{$hitId};
    die "Unknown locusId $hitId in Cofitness.hitId" if !defined $hit;
    print join("\t", $row->{orgId},
               $locusId, $gene->{sysName}, $gene->{gene}, $gene->{desc},
               $hitId, $hit->{sysName}, $hit->{gene}, $hit->{desc},
               $row->{rank}, $row->{cofit},
               exists $cons{$locusId}{$hitId} ? "TRUE" : "FALSE"
        )."\n";
}
