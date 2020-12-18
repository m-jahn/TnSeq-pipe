#!/usr/bin/perl -w

#######################################################
## keggmaplist.cgi
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Show list of KEGG maps
#
# optional: orgId, expName (multiple), which are passed on to keggmap.cgi

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;
my @expName = $cgi->param('expName');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown orgId $orgId" if $orgId && !exists $orginfo->{$orgId};

my $title = "Metabolic Maps from KEGG";
my $start = Utils::start_page($title);

print
    header,
    $start, '<div id="ntcontent">',
    h2($title),
    p("All of the metabolic maps are based on the last public release of the",
      a( {href => "http://www.genome.jp/kegg/"}, "Kyoto Encyclopedia of Genes and Genomes"),
      "(KEGG), from 2011. Gene-enzyme associations are based on the public release of KEGG and also on",
      a( {href => "http://www.jcvi.org/cgi-bin/tigrfams/index.cgi"}, "TIGRFAMS"),
      "and",
      a( {href => "http://theseed.org"}, "SEED") . "."),
    "<UL>";

my $URLsuffix = "";
$URLsuffix .= "&orgId=$orgId" if $orgId;
foreach my $expName (@expName) {
    $URLsuffix .= "&expName=$expName";
}

my $maps = $dbh->selectall_arrayref("SELECT mapId,title FROM KEGGMap ORDER BY title");
foreach my $row (@$maps) {
    my ($mapId,$mapdesc) = @$row;
    print li(a({ href => "keggmap.cgi?mapId=$mapId" . $URLsuffix },
               $mapdesc));
}

print "</UL>";
$dbh->disconnect();
Utils::endHtml($cgi);
    
