#!/usr/bin/perl -w

#######################################################
## orgAll.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page shows a table of organisms, along with how many conditions
# of fitness data are available.
#
# CGI parameters -- none

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;


my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

# main table
# gather data and slice it into an array of hashes
# my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup, COUNT(DISTINCT condition_1) AS nCond, COUNT(*) as nExp FROM Experiment GROUP BY orgId,expGroup ORDER BY orgId; },
    # { Slice => {} });
# my $cond = $dbh->selectall_arrayref(qq{SELECT *, COUNT(DISTINCT expDesc) AS nCond, COUNT(*) as nExp FROM Experiment WHERE expGroup NOT IN ('carbon source', 'nitrogen source', 'stress') GROUP BY orgId ORDER BY orgId; },
#     { Slice => {} });

my $all = $dbh->selectall_arrayref(qq{SELECT e.orgId, o.division
	FROM "Experiment" as e 
	JOIN "Organism" as o on e.orgId = o.orgId
	GROUP BY e.orgId ORDER BY o.genus, o.species, o.strain; },
    { Slice => {} });

my $count = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT expDesc) AS nCond
	FROM "Experiment"
	WHERE expGroup NOT IN ('carbon source', 'nitrogen source', 'stress') 
	GROUP BY orgId  },
    'orgId');

# my $total = $dbh->selectall_hashref(qq{SELECT *, COUNT(*) as nExp FROM Experiment GROUP BY orgId ORDER BY genus, species, strain; },
#     'orgId');

my $total = $dbh->selectall_hashref(qq{SELECT e.*, COUNT(*) as nExp FROM 
	"Experiment" as e 
	JOIN "Organism" as o on e.orgId = o.orgId
	GROUP BY e.orgId ORDER BY o.genus, o.species, o.strain; },
    'orgId');
# my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup, COUNT(DISTINCT condition_1) AS nCond, COUNT(*) as nExp FROM Experiment GROUP BY orgId ORDER BY orgId, expGroup; },
#     { Slice => {} });

my $nNitrogen = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='nitrogen source' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');
my $nCarbon = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='carbon source' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');
my $nStress = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='stress' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');



# write the title
my $title = scalar(keys %$orginfo) . " Organisms";
my $start = Utils::start_page("$title");

print
    header,
    $start, '<div id="ntcontent">',
    h2("$title");

#exit if no results
# Utils::fail($cgi, "No experiments for this organism.") if @$cond == 0;

my $nreannos = $dbh->selectall_arrayref("SELECT orgId, COUNT(*) FROM Reannotation GROUP BY orgId");
my %nreanno = map { $_->[0] => $_->[1] } @$nreannos;

#create table
my @trows = Tr({-valign => "bottom", -align => 'center'},
               th([ 'Organism', 'Division', 'Experiments']),
               th({-style => "background-color: white;"},
                  [ small('# Carbon Sources'),
                    small('# Nitrogen Sources'),
                    small('# Stress Compounds'),
                    small('# Other Conditions') ]),
               th('Updated Annotations') );
foreach my $row (@$all) {
  my $orgId = $row->{orgId};
  my $nreannoString = "&nbsp;";
  $nreannoString = a( { -href => "orgReanno.cgi?orgId=$orgId" }, $nreanno{$orgId} )
    if $nreanno{$orgId};
  push @trows, Tr({ -valign => 'top', -align => 'right' },
                  td({-align => 'left'},
                      [ a({-href=>"org.cgi?orgId=$orgId"}, $orginfo-> {$orgId}{genome}),
                        small($row->{division}) ]),
                  td({-align => 'right'},
                     a( { href => "org.cgi?orgId=$orgId"},
                          $total->{$orgId}{nExp} )),
                  td({-align => 'right', -style => 'background-color: white;'},
                     [ a({-href=>"exps.cgi?orgId=$orgId&expGroup=carbon%20source"},
                         $nCarbon->{$orgId}{nCon} || ""),
                       a({href=>"exps.cgi?orgId=$orgId&expGroup=nitrogen%20source"},
                         $nNitrogen->{$orgId}{nCon} || ""),
                       a({href=>"exps.cgi?orgId=$orgId&expGroup=stress"},
                         $nStress->{$orgId}{nCon} || ""),
                       $count->{$orgId}{nCond} || "" ]),
                  td({-align => 'right'}, $nreannoString)
		 );
}

print table({cellspacing => 0, cellpadding => 3}, @trows);
my $nTotExps = 0;
foreach my $row (@$all) {
    $nTotExps += $total->{$row->{orgId}}->{nExp};
}
print p("Total number of genome-wide fitness experiments: ", Utils::commify($nTotExps));

# display number of genes that have data out of total genes
# gather number of genes and data
# my $numGenes = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM Gene;}, undef);
# my $numData = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM GeneFitness;}, undef);
# print $cgi->p("Fitness data for $numData genes of $numGenes genes.");


#display taxonomic information and link
# if ((defined $orginfo->{$orgId}{taxonomyId}) && ($orginfo->{$orgId}{taxonomyId} ne "")) {
	# print $cgi->p($cgi->a({href => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$orginfo->{$orgId}{taxonomyId}"}, "NCBI Taxonomy"));
# }

#file of experimental data - generated by createExpData.cgi
# my $data = `./createExpData.cgi orgId=$orgId`;
# print $cgi->p($cgi->a({href => "download.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a minute or so to load once clicked.");
# print $cgi->p($cgi->a({href => "createExpData.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a few seconds to load once clicked.");

$dbh->disconnect();
Utils::endHtml($cgi);
