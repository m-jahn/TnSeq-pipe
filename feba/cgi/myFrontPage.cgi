#!/usr/bin/perl -w
#######################################################
## myFrontPage.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu), and 
## Morgan Price
#######################################################
#
# CGI parameters: none

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
# my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my @orgOptions = ("");
my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
my %orgLabels = ("" => join(" ", "All", scalar(@orginfo), "organisms"));
foreach my $hash (@orginfo) {
    my $orgId = $hash->{orgId};
    push @orgOptions, $orgId;
    $orgLabels{$orgId} = $hash->{genome};
}
my $start = Utils::start_page("Fitness Browser");

print header, $start,
 #    start_html(
 #      $start,
	# -title =>"Fitness Browser",
	# # -style => {-code => $style},
	# -author=>'wjshaoATberkeley.edu, victorialoATberkeley.edu',
	# -meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),

    
div({-id=>"ntcontent"},
  h2("Fitness Browser"),
  Utils::site_intro_text(),
  h3("Search"),
  table({class=>"frontpage"}, 
    $cgi->Tr(
      td({class=>"searchbar"},
        a({href=>"geneSearch.cgi"},'by Gene<BR>',img({src=>"../images/genesearch.png"}))
        ),
      td({class=>"searchbar"},
        a({href=>"blastSearch.cgi"},'by Sequence<BR>',img({src=>"../images/bysequence.png"}))
        ),
      td({class=>"searchbar"},
        a({href=>"expSearch.cgi"},'by Condition<BR>',img({src=>"../images/bycondition.png"}))
        ),
      td({class=>"searchbar"},
        a({href=>"orgAll.cgi"},'Organisms<BR>',img({src=>"../images/organisms.png"}))
        ),
      td({class=>"searchbar"},
        'See Conditions:<BR><BR>', a({href=>"cond.cgi?expGroup=carbon%20source", title=>"Carbon Sources"},"Carbon Sources<BR>"), a({href=>"cond.cgi?expGroup=nitrogen source", title=>"Nitrogen Sources"},"Nitrogen Sources<BR>"), a({href=>"cond.cgi?expGroup=stress", title=>"Stresses"},"Stresses"),
        )
      )
    ),
    h3("Examples"),
    table({class=>"frontpage"},
      Tr(
        th("Genes"),
        td(a({href=>"myFitShow.cgi?orgId=Keio&gene=18086&help=1"},"Top Phenotypes<BR>", img({src=>"../images/heatmap1.png"}))),
        td(a({href=>"cofit.cgi?orgId=Keio&locusId=18086&help=1"},"Top Cofitness<BR>", img({src=>"../images/cofit.png"}))),
        td(a({href=>"genesFit.cgi?orgId=Keio&showAll=0&locusId=18086&locusId=14904&addgene=uvrC&help=1"},"Compare Genes<BR>", img({src=>"../images/heatmapN.png"})))

      ), Tr(
        th("Experiments"),
        td(a({href=>"exp.cgi?orgId=Keio&expName=set2IT043&show=specific&help=1"},"Specific Phenotypes<BR>", img({src=>"../images/cofit.png"}))),
        td(a({href=>"compareExps.cgi?orgId=Keio&expName2=set1IT031&expName1=set1IT003&Go=Go&help=1"},"Compare Experiments<BR>", img({src=>"../images/CmpExps.png"}))),
        td(a({href=>"compareExps.cgi?orgId=Keio&expName2=set1IT031&expName1=set1IT003&outlier=lowy&minabs=1&submit=Go&help=1"},"Outliers<BR>", img({src=>"../images/outliers.png"})))
      ), Tr(
        th("Multi-Organism"),
        td(a({href=>"orthCond.cgi?expGroup=stress&condition1=Cisplatin&help=1"},"Overview of a Condition<BR>", img({src=>"../images/cond_overview.png"}))),
        td(a({href=>"orthFit.cgi?orgId=Keio&locusId=14904&expGroup=stress&condition1=Cisplatin&help=1"},"Conserved Phenotypes<BR>", img({src=>"../images/CmpFit.png"}))),
        td(a({href=>"cofitCons.cgi?orgId=Keio&locusId=14904&hitId=18086&help=1"},"Conserved Cofitness<BR>", img({src=>"../images/CmpCofit.png"})))
      )
    ),
  h6(q(Developed by <A HREF="http://morgannprice.org/">Morgan Price</A>, Victoria Lo, and Wenjun Shao
         in the <A HREF="http://genomics.lbl.gov">Arkin lab</A>.
         Please report any bugs to <A HREF="mailto:funwithwords26@gmail.com">Morgan</A>.)),
  );

      
 #    # Gene search form
 #    div({-style => "float:left; vertical-align: top;"},
	# h3("Search by Gene Name"),
	# start_form( -name    => 'input', -method  => 'GET', -action  => 'myFitShow.cgi' ),
	# # drop down list of species
	# p("Choose organism:",
	#   popup_menu( -name    => 'orgId', -values  => \@orgOptions, -labels  => \%orgLabels, -default => $orgOptions[0])),
	# p("Enter gene name:",
	#   textfield( -name      => 'gene', -size => 20, -maxlength => 100 ),
	#   br(),
	#   small(qq{examples: "Shewana3_0001" or "recA"})),
	# p(submit("Find gene")),
	# end_form ),
 #    # Experiment search form
 #    div({-style => "padding-left: 10px; padding-right: 10px; float: right; vertical-align: top; background-color: rgb(240,240,240);"},
	# h3("Find Experiments"),
	# start_form(-name => 'byExp', -method => 'GET', -action => 'exps.cgi'),
	# p("Choose organism:",
	#   popup_menu(-name => 'orgId', -values => \@orgOptions, -labels => \%orgLabels, -default => $orgOptions[0])),
	# p("Enter condition:", textfield(-name => 'query', -size => 20, -maxlength => 100), submit("Go"),
	#   br(),
	#   small(qq{examples: "cisplatin" or "glucose"})),
	# p("Or show",submit("All experiments"),"for one organism"),
	# p("Or see all conditions:", br(),
	#   a({ href => "cond.cgi?expGroup=carbon source" }, "carbon sources," ),
	#   a({ href => "cond.cgi?expGroup=nitrogen source" }, "nitrogen sources," ),
	#   "or",
	#   a({ href => "cond.cgi?expGroup=stress" }, "stresses" )),
 #        end_form),
 #    div({-style => "clear: right"}),
 #    h3(qq(Search by Gene Sequence)),
 #    start_form( -name    => 'input', -method  => 'POST', -action  => 'mySeqSearch.cgi' ),
 #    p("Choose query type: ", 
 #      popup_menu( -name => 'qtype', -values => [ 'protein', 'nucleotide' ], -default => 'protein' )),
 #    p("Enter query sequence: ",
 #      br(),
 #      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
 #    p("How many hits to show: ",
 #      popup_menu( -name => 'numHit', -values  => [ 5,10,25,50,100 ], -default => 25 )),
 #    p(submit("Search"),
 #      reset() ),
 #    end_form,
    
    end_html;

$dbh->disconnect();

# END
