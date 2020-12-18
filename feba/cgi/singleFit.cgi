#!/usr/bin/perl -w
#######################################################
## singleFit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu), and 
## Morgan Price
#######################################################
#
# Redirected to from myFitShow if there is only one match.
# Created for reorganization and to fix issues with tabs.
# 
# Required CGI parameters:
# orgId -- which organism to search in
# locusId -- which locus to search on
#
# Optional CGI parameters: 
# showAll -- 1 if showing all fitness values instead of just the most extreme ones
# help -- 1 if on help/tutorial mode
 
use strict;

use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use HTML::Entities;

use lib "../lib";
use URI::Escape;
use Utils;

my $cgi=CGI->new;
print $cgi->header;
my $style = Utils::get_style();

my $orgSpec = $cgi->param('orgId') || "";
my $locusId = $cgi->param('locusId');
my $showAll = $cgi->param('showAll') ? 1 : 0;
my $help = $cgi->param('help') || "";
my $locusSafe = HTML::Entities::encode($locusId);

# connect to database
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId argument" if $orgSpec ne "" && !exists $orginfo->{$orgSpec};

my $start = Utils::start_page("Fitness for $locusSafe ($orgSpec)");

my $query = qq{SELECT orgId, locusId, sysName, desc, gene, type FROM Gene
		WHERE locusId = ? };
my $hits;
if ($orgSpec) {
  $query .= " AND orgId = ?";
  $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $locusId, $orgSpec);
} else {
  $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $locusId);
}

# and add another column for whether there is fitness data
foreach my $gene (@$hits) {
  $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
}

if (@$hits == 0) {
  print $start,'<div id="ntcontent">',
    $cgi->h3("No gene found for",  $locusSafe,
             (exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : ""));

} else {
  # just 1 hit
  my $gene = $hits->[0];
  my $orgId = $gene->{orgId};
  my $locusId = $gene->{locusId};

  my $idShow = $gene->{sysName} || $gene->{locusId};
  my $start = Utils::start_page("Fitness data for $idShow in $orginfo->{$orgId}{genome}");
  my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,$gene->{type},"fitness");
  print
    $start,
    $tabs,
    h2("Fitness data for $idShow in " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
    p(Utils::gene_link($dbh, $gene, "lines"));

  if ($hits->[0]{has_fitness} == 0) {
    print p("Sorry, no fitness data for $idShow.",
              "This gene might not have mapped insertions due to chance,",
              "or its mutants might be at low abundance after recovery from the freezer,",
              "or its mutants might grow poorly in the conditions that were used to make the mutant library.");
    my $exps = $dbh->selectall_arrayref("SELECT * from Experiment WHERE orgId = ? LIMIT 1", { Slice => {} }, $orgId);
    if (scalar(@$exps) > 0) {
      my $expName = $exps->[0]{expName};
      print p("You can see the insertions near $idShow that have per-strain fitness values",
              a({href => "strainTable.cgi?orgId=$orgId&locusId=$locusId&expName=$expName"}, "here") . ".",
              "Insertions that are at low abundance after recovery from the freezer",
              "or that lack unique barcodes are not shown.");
    }
  } else {
      # show fitness data for gene
      my @fit = @{ $dbh->selectall_arrayref("SELECT expName,fit,t from GeneFitness where orgId=? AND locusId=?",
                                            { Slice => {} },
                                            $orgId, $locusId) };
      my $nTotValues = scalar(@fit);
      die "Unreachable" if $nTotValues == 0;
      my $limitRows = $showAll ? $nTotValues : 20;
      my $minAbsFit;
      if ($nTotValues > $limitRows) {
        # subset the rows
        @fit = sort { abs($b->{fit}) <=> abs($a->{fit}) } @fit;
        @fit = @fit[0..($limitRows-1)];
        $minAbsFit = abs($fit[$#fit]{fit});
      }

      # and get metadata about experiments
      my $expinfo = Utils::expinfo($dbh,$orgId);

      if ($showAll) {
        @fit = sort { Utils::CompareExperiments($expinfo->{$a->{expName}}, $expinfo->{$b->{expName}}) } @fit;
      } else {
        @fit = sort { $a->{fit} <=> $b->{fit} } @fit;
      }

      if ($help) {
          print qq[<div class="helpbox">
		<b><u>About this page:</u></b><BR><ul>
		<li>View the experiments that have the strongest phenotypes for this gene. (Or choose to view all phenotypes with the link below.)</li>
		<li>To get to this page, search for any gene and click on the "Fitness" tab.</li>
		<li>Use the "Add gene box" above to compare this gene's fitness to another one (try "uvrB").</li>
		<li>Hover over blue links for more information.</li>
                <li>For more about fitness values, see the <A HREF="help.cgi#fitness">help page</A>.
		</ul></div>];
	}

	if ($showAll) {
          print $cgi->p("All " . scalar(@fit) . " fitness values, sorted by group and condition");
	} else {
          if (defined $minAbsFit) {
            $minAbsFit = sprintf("%.1f", $minAbsFit);
            print $cgi->p("Top $limitRows experiments with the strongest phenotypes (|fitness| &ge; $minAbsFit)");
          } else {
            print $cgi->p("All " . scalar(@fit) . " fitness values, sorted by value");
          }
	}

      print
	    qq[<div style="position: relative;"><div class="floatbox">],
              start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
                "<br>Add gene: ",
                  hidden( -name => 'orgId', -value => $orgId, -override => 1 ),
                    hidden( -name => 'showAll', -value => $showAll, -override => 1  ),
                      hidden( -name => 'locusId', -value => $locusId, -override => 1 ),
                        textfield( -name => 'addgene', -default => "", -override => 1, -size => 20, -maxLength => 100 ),
                          # submit('Go'),
                          "<button type='submit'>Go</button>",
                            end_form,
                              qq[</P></div></div>];


      my $option = "Or see ";
      if ($showAll == 0) {
        my $showAllDest = qq(singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=1);
        $option .= qq[<a href="$showAllDest">all $nTotValues experiments</a>];
      } else {
        my $showFewDest = qq(singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0);
        $option .= qq[<a href="$showFewDest">strongest phenotypes</a>];
      }
      print p($option, " or ",
              a({ -href => "heatmap.cgi?orgId=$orgId&r=$locusId" }, "choose conditions"));
      my @out = (); # specifiers for HTML rows for each fitness value
      my $lastGroup = ""; # only enter the first time
      foreach my $fitrow (@fit) {
        my $expName = $fitrow->{expName};
        my $exp = $expinfo->{$expName};
        my $group = $exp->{expGroup};
        my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$locusId&expName=$expName";
        $strainUrl .= "&help=1" if $help;
        push @out, join(" ",
                        td($group eq $lastGroup ? "" : $group),
                        td(a({ 
                              title => "$expName: $exp->{expDescLong}",
                              href => "exp.cgi?orgId=$orgId&expName=$expName" },
                             $exp->{expDesc})),
                        td( { -bgcolor => Utils::fitcolor($fitrow->{fit}) },
                            a({ -href => $strainUrl,
                                -title => "per-strain data",
                                -style => "color:rgb(0,0,0)" },
                              sprintf("%.1f", $fitrow->{fit}) ) ),
                        td( sprintf("%.1f", $fitrow->{t}) ),
                        td(a({ title => "Compare to data from similar experiments or orthologs",
                               href => "orthFit.cgi?orgId=$orgId&locusId=$locusId"
                               . "&expGroup=" . uri_escape($exp->{expGroup})
                               . "&condition1=" . uri_escape($exp->{condition_1}) }),
                           "compare") );
        $lastGroup = $group if $showAll;
      }
      my $relsize = $showAll ? "70%" : "100%";
      print $cgi->table( { cellspacing => 0, cellpadding => 3, },
                         $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                                  $cgi->th( [ 'group', 'condition','fitness','t score', '&nbsp;' ] ) ),
                         $cgi->Tr({-align=>'left',-valign=>'top',-style=>"font-size: $relsize"}, \@out ) );

    print "<br><br>";
    } #  end if just 1 hit
}

$dbh->disconnect();
Utils::endHtml($cgi);
