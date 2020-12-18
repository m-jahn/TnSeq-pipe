#!/usr/bin/perl -w
#######################################################
## orthFit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Given Group, Condition_1, and an anchor gene, list fitness values
# of the gene and its orthologs for relevant conditions.
# This is useful even if there are no orthologs because it
# shows other replicates or other concentrations.
#
# Required CGI parameters:
# orgId, locusId -- these specify the anchor gene
# expGroup and condition1 -- these specify the condition.
#	Note that condition1 may be empty but it must still be present
#	(and is used to choose which to show, i.e. matchingon empty condition1).
# Optional: help -- 1 if on help/tutorial mode

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "no orgId parameter";
my $locusId = $cgi->param('locusId') || die "no locusId parameter";
# Allow expGroup missing because occasionally an experiment does not have an expGroup
# (Should that be an error during loading?)
my $expGroup = $cgi->param('expGroup') || "";
my $condition1 = $cgi->param('condition1');
die "no condition1 parameter" unless defined $condition1;
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId" unless exists $orginfo->{$orgId};
my $gene = $dbh->selectrow_hashref("SELECT * from Gene where orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
die "No such gene" unless defined $gene->{locusId};

my $orths = Utils::get_orths($dbh, $orgId, $locusId);
my @orths = sort { $b->{ratio} <=> $a->{ratio} } values(%$orths);
my @genes = ($gene);
push @genes, @orths;
# ensure that condition1 is valid
my ($cnt) = $dbh->selectrow_array(qq{SELECT count(*) from Experiment WHERE orgId = ? AND expGroup = ? AND condition_1 = ?},
				  {}, $orgId, $expGroup, $condition1);
die "No experiments for specified expGroup and condition1\n" unless $cnt > 0;

my $showId = $gene->{sysName} || $gene->{locusId};
my $title = "$showId and its orthologs: fitness in $expGroup $condition1";
my $start = Utils::start_page("$title");

print header,
	$start, '<div id="ntcontent">',
    h2($title);

if ($help) {
        print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View all of the fitness data for a gene and its orthologs in a given condition. </li>
        <li>To get to this page, search for a gene, visit the "Fitness" tab, and click on the "compare" link for a condition of interest.</li> 
        <li>Or, search for an experiment, visit the "Specific" tab, and click on the "compare" link for a gene of interest.</li>
        <li>Click on links (including the fitness values) to view more.</li>
        </ul></div><BR><BR>];
    }


    print div({-style => "clear: right"});

my @headings = (a({title=>"score ratio for ortholog: blast_score_of_alignment / self_score"},'Ratio'),
		 qw{Organism Gene Name Description Experiment Fitness t} );
my @trows = ( $cgi->Tr({ -valign=> 'top', -align => 'center'},
		       map { th($_) } @headings) );
my $nSkipOrth = 0;

my $shade = 0;
foreach my $o (@genes) {
  my $data = $dbh->selectall_arrayref(qq{SELECT * from Experiment JOIN GeneFitness USING  (orgId,expName)
                                           WHERE orgId=? AND locusId=? AND expGroup=? AND condition_1=? ORDER BY fit},
                                      { Slice => {} }, $o->{orgId}, $o->{locusId}, $expGroup, $condition1);
  if (@$data == 0 && $o->{orgId} ne $orgId) {
    $nSkipOrth++;
    next;
  }
  my $first = 1;
  $shade++;
  foreach my $row (@$data) {
    my $ratio = $o->{orgId} eq $orgId ? "&mdash;" : sprintf("%.2f",$o->{ratio});
    my $orgShort = "";

    if ($first) {
      my $d = $orginfo->{$o->{orgId}};
      my $short = $d->{genome}; $short =~ s/^(.)[^ ]+/$1./;
      $orgShort = a({href => "org.cgi?orgId=$o->{orgId}", title => $d->{genome}}, $short);
    }
    my $strainUrl = "strainTable.cgi?orgId=$o->{orgId}&locusId=$o->{locusId}&expName=$row->{expName}";
    $strainUrl .= "&help=1" if $help;
    push @trows,
      $cgi->Tr( { -valign => 'top', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF" },
                td($first ?  $ratio : ""),
                td($orgShort),
                td($first ? Utils::gene_link($dbh, $o, "name", "myFitShow.cgi") : ""),
                td($first ? $o->{gene} : ""),
                td($first ? Utils::gene_link($dbh, $o, "desc", "domains.cgi") : ""),
                td(a({href => "exp.cgi?orgId=$o->{orgId}&expName=$row->{expName}"},
                     $row->{expDesc})),
                td({ -bgcolor => Utils::fitcolor($row->{fit}) },
                   a({ -href => $strainUrl,
                       -title => "per-strain data",
                       -style => "color:rgb(0,0,0)" },
                     sprintf("%.1f", $row->{fit}) ) ),
                td(sprintf("%.1f", $row->{t})) );
    $first = 0;
  }
}
print table({ cellpadding => 3, cellspacing =>0}, @trows);
print p("$nSkipOrth orthologs are not shown because they lack fitness data for this condition (or they lack data entirely)")
  if $nSkipOrth;
print p("$showId is not a protein-coding gene, so it has no orthologs.") if $gene->{type} != 1;
print p(a({href => "orthCond.cgi?expGroup=" . uri_escape($expGroup) . "&condition1=" . uri_escape($condition1) },
      "Specific phenotypes for $expGroup $condition1 across organisms"));
print p(a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Show all homologs"))
  if $gene->{type} == 1;
print p(a({-href => "cmpbrowser.cgi?anchorOrg=$orgId&anchorLoci=$locusId"}, "Comparative fitness browser"));

$dbh->disconnect();
Utils::endHtml($cgi);
