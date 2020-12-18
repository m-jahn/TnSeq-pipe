#!/usr/bin/perl -w
#######################################################
## pathwaysOrg.cgi -- list MetaCyc pathways for an organism
##
## Copyright (c) 2018 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

print
  header,
  Utils::start_page("MetaCyc Pathways for $orginfo->{$orgId}{genome}"),
  h2("Metacyc Pathways for",
     a({ -href => "org.cgi?orgId=$orgId" }, $orginfo->{$orgId}{genome}));

my $paths = $dbh->selectall_arrayref(qq{ SELECT * FROM MetacycPathway
                                         JOIN MetacycPathwayCoverage USING (pathwayId)
                                         WHERE orgId = ? },
                                     { Slice => {} }, $orgId);
foreach my $path (@$paths) {
  $path->{score} = Utils::MetacycPathwayScore($path->{nSteps}, $path->{nFound});
}
my @sorted = sort { $b->{score} <=> $a->{score}
                    || $a->{nSteps} <=> $b->{nSteps}
                    || $a->{pathwayName} cmp $b->{pathwayName} } @$paths;

my @trows = ();
push @trows, Tr(th("Pathway"), th("Steps Found"));

foreach my $path (@sorted) {
  my $pathId = $path->{pathwayId};
  push @trows, Tr( td( a( { -href => "pathway.cgi?orgId=$orgId&pathwayId=$pathId" },
                          $path->{pathwayName } ) ),
                   td( $path->{nFound}, "/", $path->{nSteps} ) );
}
print
  table( { cellspacing => 0, cellpadding => 3 }, @trows),
  p(small("Only pathways with at least one candidate gene are shown"));

$dbh->disconnect();
Utils::endHtml($cgi);
