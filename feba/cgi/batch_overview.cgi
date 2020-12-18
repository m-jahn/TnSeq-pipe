#!/usr/bin/perl -w
#######################################################
## batch_overview.cgi -- show best hits for a gene from a batch comparison
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# jobId -- the job identifier

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use Batch;

my $cgi=CGI->new;

my $jobId = $cgi->param('jobId');
die "Much specify jobId" unless $jobId ne "";

my $bdb = Batch::get_batch_dbh($jobId);
my $jobname = Batch::get_job_name($jobId,$bdb);

my $stats = $bdb->selectrow_hashref(qq{SELECT COUNT(*) AS nProt,
                                       SUM(hasSpecific OR hasCofit) AS nLink,
                                       SUM(isHypo) AS nHypo,
                                       SUM(isHypo AND (hasSpecific OR hasCofit)) AS nHypoLink
                                       FROM Query});

my ($nWithBestHit) = $bdb->selectrow_array(qq{SELECT count(DISTINCT queryId) FROM BestHit});
my ($nWithBBH) = $bdb->selectrow_array(qq{SELECT count(DISTINCT queryId) FROM BestHit WHERE isBBH});

print header,
    Utils::start_page("Fitness BLAST - $jobname"),
    '<div id="ntcontent">',
    h2("Fitness BLAST for $jobname"),
    h3("Statistics"),
    ul(li("Best hits for $nWithBestHit of $stats->{nProt} proteins"),
       ul(li("$nWithBBH with orthologs")),
       li("Functional links for $stats->{nLink} of $stats->{nProt} proteins"),
       ul(li("and for $stats->{nHypoLink} of $stats->{nHypo} hypotheticals"))),
    h3("Browse Proteins"),
    start_form( -name => 'input', -method => 'GET', -action => 'batch_genes.cgi' ),
    hidden( -name => 'jobId', -value => $jobId, -override => 1 ),
    p("Types of proteins:", radio_group(-name => 'hypo',
                               -values => ["All","Hypotheticals"],
                               -default => "All")),
    p("With functional link: ", radio_group(-name => 'link',
                                            -values => ["Specific Phenotype","Cofit","Either","Not needed"],
                                            -default => "Either")),
    p(submit( { -style => 'position: absolute; text-align: left; float: none; display:inline;' , -name => "Browse" })),
    end_form,
    p(" &nbsp; "),
    p(" &nbsp; "),
    h3("Search by Identifier or Keyword"),
    start_form( -name => 'input', -method => 'GET', -action => 'batch_genes.cgi' ),
    hidden( -name => 'jobId', -value => $jobId, -override => 1 ),
    p(textfield(-name => "search", -size => 30, -maxlength => 100),
      "&nbsp;",
      submit( { -style => 'position: absolute; text-align: left; float: none; display:inline;' , -name => "Search" })),
    end_form;

my $dir = "../job_data/$jobId"; # get_batch_dbh ensures its validity

print small(br(),
            p("The analysis was run on ",
              `date -r ../job_data/$jobId/hits`,
              "- you can",
              a({-href=>"batch_blast.cgi?jobId=$jobId"}, "rerun the analysis")));
Utils::endHtml($cgi);

