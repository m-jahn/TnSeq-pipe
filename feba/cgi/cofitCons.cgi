#!/usr/bin/perl -w
#######################################################
## cofit.cgi -- examine cofitness of orthologs of two genes
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId, locusId, hitId
# Optional: help -- 1 if on help/tutorial mode

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

sub cofitStrings($$$$); # dbh, orgId, locus1, locus2 => cofitness, rank
my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no orgId";
my $locusId = $cgi->param('locusId') || die "no locusId";
my $hitId = $cgi->param('hitId') || die "no hitId";
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId" unless exists $orginfo->{$orgId};

my $gene1 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?;",
				    {}, $orgId, $locusId);
die "Invalid locusId" unless defined $gene1->{locusId};
my $gene2 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?;",
				    {}, $orgId, $hitId);
die "Invalid hitId" unless defined $gene2->{locusId};

my $show1 = $gene1->{sysName} || $gene1->{locusId};
my $show2 = $gene2->{sysName} || $gene2->{locusId};

my $orth1 = Utils::get_orths($dbh,$orgId,$locusId);
my $orth2 = Utils::get_orths($dbh,$orgId,$hitId);
my @orth1 = sort { $b->{ratio} <=> $a->{ratio} } values(%$orth1);
my @orth1b = grep { exists $orth2->{$_->{orgId}} } @orth1;

my @trows = ( Tr({-valign => 'top', -align => 'center'},
		 th(["Organism",
		     a({title=>"score ratio for ortholog 1: blast_score_of_alignment / self_score"},'Ratio1'),
		     "Gene1", "Name1", "Description1",
		     a({title=>"score ratio for ortholog 2: blast_score_of_alignment / self_score"},'Ratio2'),
		     "Gene2", "Name2", "Description2", "Cofit", "Rank" ])),
	     Tr({-valign => 'top', -align => 'left'},
		 td([ $cgi->a({ -href => "org.cgi?orgId=". $orginfo->{$orgId}{orgId},
                                -title => $orginfo->{$orgId}{division} },
                              $orginfo->{$orgId}{genome} ),
		      "1.0",
		      Utils::gene_link($dbh, $gene1, "name", "myFitShow.cgi"),
		      $gene1->{gene},
		      Utils::gene_link($dbh, $gene1, "desc", "domains.cgi"),

		      "1.0",
		      Utils::gene_link($dbh, $gene2, "name", "myFitShow.cgi"),
		      $gene2->{gene},
		      Utils::gene_link($dbh, $gene2, "desc", "domains.cgi"),

		      &cofitStrings($dbh,$orgId,$locusId,$hitId) ])) );

my $nBoth = scalar(@orth1b);
my $n1only = scalar(@orth1) - $nBoth;
my $n2only = scalar(keys %$orth2) - $nBoth;

foreach my $o1 (@orth1b) {
    my $orthOrg = $o1->{orgId};
    my $o2 = $orth2->{$orthOrg};
    push @trows,
      Tr( {-valign => 'top', -align => 'left'},
          td([ $cgi->a({ -href => "org.cgi?orgId=". $orginfo->{$orthOrg}{orgId},
                         -title => $orginfo->{$orthOrg}{division} },
                       $orginfo->{$orthOrg}{genome} ),
               sprintf("%.2f", $o1->{ratio}),
               Utils::gene_link($dbh, $o1, "name", "myFitShow.cgi"),
               $o1->{gene},
               Utils::gene_link($dbh, $o1, "desc", "domains.cgi"),
               sprintf("%.2f", $o2->{ratio}),
               Utils::gene_link($dbh, $o2, "name", "myFitShow.cgi"),
               $o2->{gene},
               Utils::gene_link($dbh, $o2, "desc", "domains.cgi"),
               &cofitStrings($dbh, $orthOrg, $o1->{locusId}, $o2->{locusId}) ]));
}

my $title = "Conservation of cofitness between $show1 and $show2 in $orginfo->{$orgId}{genome}";
my $start = Utils::start_page("$title");
print
	$start, '<div id="ntcontent">',
    h2($title);

if ($help) {
        print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>See if two genes have similar fitness patterns (high <A HREF="help.cgi#cofitness">cofitness</A>) across organisms. </li>
        <li>To get to this page, search for a gene, visit the "Cofit" tab, and click on "compare".</li>
        <li>To see the actual fitness patterns in one organism, click on the cofitness value. On that page, you can request a scatterplot.</li>
        </ul></div>];
    }

    print h3("$nBoth genomes with putative orthologs of both genes"),
    table({-cellspacing => 0, cellpadding => 3}, @trows);
print p("Not shown: $n1only genomes with orthologs for $show1 only; $n2only genomes with orthologs for $show2 only");

# then print table of the two genes
# then print graceful error message if either lacks orthologs
# then print graceful error message if no organisms with orthologs for both (@orth1b empty)
# then for each in 1 with orth1b, look for cofitness...

Utils::endHtml($cgi);

sub cofitStrings($$$$) {
    my ($dbh,$orgId,$locus1,$locus2) = @_;
    my ($cofit,$rank) = $dbh->selectrow_array("SELECT cofit,rank FROM Cofit where orgId=? AND locusId=? AND hitId=?",
					      {}, $orgId, $locus1, $locus2);
    my $URL = "genesFit.cgi?orgId=$orgId&locusId=$locus1&locusId=$locus2";
    return (a({href => $URL}, sprintf("%.2f", $cofit)), $rank) if defined $cofit;
    # else !defined cofit
    # get min cofit and max rank
    my $hasData1 = Utils::gene_has_fitness($dbh,$orgId,$locus1);
    my $hasData2 = Utils::gene_has_fitness($dbh,$orgId,$locus2);
    if ($hasData1 && $hasData2 ) {
	($cofit,$rank) = $dbh->selectrow_array("SELECT min(cofit), max(rank) FROM Cofit where orgId=? AND locusId=? GROUP BY orgId,locusId;",
					       {}, $orgId, $locus1);
	return (a({href => $URL, title => sprintf("under %.2f",$cofit)}, "low"), "&gt; $rank") if defined $cofit;
	# else
	return (a({title => "no cofitness for this organism"}, "&mdash;"), a({title => "no cofitness for this organism"}, "&mdash;"));
    }
    # else
    return (a({title => "no data for either"}, "&mdash;"), a({title => "no data for either"}, "&mdash;"))
	if (!$hasData1 && !$hasData2);
    return (a({title => "no data for gene 1"}, "&mdash;"), a({title => "no data for gene 1"}, "&mdash;"))
	if !$hasData1 && $hasData2;
    return (a({title => "no data for gene 2"}, "&mdash;"), a({title => "no data for gene 2"}, "&mdash;"))
	if $hasData1 && !$hasData2;
    die "Unreachable";
}

