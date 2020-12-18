#!/usr/bin/perl -w
#######################################################
## geneSearch.cgi
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
my $start = Utils::start_page("Fitness Browser - Exp Search");

print header, $start,    
div({-id=>"ntcontent"},
    h2("Search Experiments or Conditions"),
    Utils::site_intro_text(),
    start_form(-class => 'search', -name => 'byExp', -method => 'GET', -action => 'exps.cgi'),
    p("Choose organism:",
      popup_menu(-name => 'orgId', -values => \@orgOptions, -labels => \%orgLabels, -default => $orgOptions[0])),
    p("Enter condition:", textfield(-name => 'query', -size => 20, -maxlength => 100), br(), small(qq{examples: "cisplatin" or "glucose"}), submit("Go"),
    ),
    p("Or see all conditions:", br(),
      a({ href => "cond.cgi?expGroup=carbon%20source" }, "carbon sources," ),
      a({ href => "cond.cgi?expGroup=nitrogen%20source" }, "nitrogen sources," ),
      "or",
      a({ href => "cond.cgi?expGroup=stress" }, "stresses" )),
    end_form,

    h6(q(Developed by Victoria Lo, Wenjun Shao, and
         <A HREF="http://morgannprice.org/">Morgan Price</A>.
         Please report any bugs to <A HREF="mailto:funwithwords26@gmail.com">Morgan</A>.)),
    ),
    end_html;

$dbh->disconnect();

# END
