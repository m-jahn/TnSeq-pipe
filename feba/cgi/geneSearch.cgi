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
my $start = Utils::start_page("Fitness Browser - Gene Search");

print header, $start,
    
div({-id=>"ntcontent"},
    h2("Search by Gene Name"),
    Utils::site_intro_text(),
    start_form( -class => "search", -name    => 'input', -method  => 'GET', -action  => 'myFitShow.cgi' ),
    # drop down list of species
    p("Choose organism:",
      popup_menu( -name    => 'orgId', -values  => \@orgOptions, -labels  => \%orgLabels, -default => $orgOptions[0])),
    p("Enter gene name or domain name:",
      textfield( -name      => 'gene', -size => 20, -maxlength => 100 ),
      br(),
      small("examples:",
            map { a({href => "myFitShow.cgi?gene=$_"}, $_) }
            qw{Shewana3_0001 recA homoserine HisG_c ec:5.5.1.2 ko:K05350} )
    ),
    p(submit("Find gene")),
    end_form,
    h6(q(Developed by Victoria Lo, Wenjun Shao, and
               <A HREF="http://morgannprice.org/">Morgan Price</A>.
               Please report any bugs to <A HREF="mailto:funwithwords26@gmail.com">Morgan</A>.))
    ),
    end_html;

$dbh->disconnect();

# END
