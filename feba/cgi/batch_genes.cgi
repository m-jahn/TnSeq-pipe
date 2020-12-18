#!/usr/bin/perl -w
#######################################################
## batch_genes.cgi -- show a list of query genes that are part of a batch comparison
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# jobId -- the job identifier
#
# Various ways to list the genes:
# search -- by query identifier or word(s) in its description
# link -- use "spec" "Specific..." or "cofit" to see only genes with that type of link, or "either".
#	"any" or "not needed" is ignored. Case is ignored.
# hypo -- show only genes that lack specific annotations. Ignored if search is set.
#	"All" means show all.
#
# optional: index -- the first index of the matches to show (defaults to 0)

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use HTML::Entities;

use lib "../lib";
use Utils;
use Batch;

my $cgi=CGI->new;
my $style = Utils::get_style();

my $jobId = $cgi->param('jobId');
my $searchSpec = $cgi->param('search');

my %validLink = map {$_ => 1} qw{either cofit spec any};
my $linkSpec = $cgi->param('link') || "any";
$linkSpec = lc($linkSpec);
$linkSpec = "spec" if $linkSpec =~ m/^specific/;
$linkSpec = "any" if $linkSpec =~ m/^not/;
die "Invalid link $linkSpec" unless exists $validLink{$linkSpec};

my $onlyHypo = $cgi->param('hypo') || 0;
$onlyHypo = 0 if $onlyHypo eq "All";
my $index = $cgi->param('index') || 0;
die "Invalid index" unless $index =~ m/^\d+$/;

my $bdb = Batch::get_batch_dbh($jobId);
my $jobname = Batch::get_job_name($jobId,$bdb);

my $maxToShow = 25;

my @where = ();
my $title;
my $linkTitle = "";
if (defined $searchSpec && $searchSpec =~ m/\S/) {
    # make the search safe to include in SQL
    $searchSpec =~ s/^\s+//;
    $searchSpec =~ s/\s+$//;
    $searchSpec =~ s/[\"\n\r]//g;
    
    push @where, qq{queryId = "$searchSpec"
                   OR desc = "$searchSpec" OR desc LIKE "$searchSpec %"
                   OR desc LIKE "% $searchSpec" OR desc LIKE "% $searchSpec %"};
    # note that other specifiers are ignored if using search
    my $searchSpecShow = HTML::Entities::encode($searchSpec);
    $title = qq{Search Results for "$searchSpecShow" in $jobname};
} else {
    if ($linkSpec eq "spec") {
        push @where, "hasSpecific=1";
        $linkTitle = "with specific phenotypes";
    } elsif ($linkSpec eq "cofit") {
        push @where, "hasCofit=1";
        $linkTitle = "with cofitness";
    } elsif ($linkSpec eq "either") {
        push @where, "(hasSpecific=1 OR hasCofit=1)";
        $linkTitle = "with a functional link";
    }
    if ($onlyHypo) {
        push @where, "isHypo=1";
        $title = qq{Hypothetical Proteins in $jobname};
    } else {
        $title = qq{Proteins in $jobname};
    }
    $title .= " with orthologs $linkTitle" if $linkTitle ne "";
}
my $whereClause = @where > 0 ? "WHERE " . join(" AND ", @where) : "";
my $qgenes = $bdb->selectall_arrayref(
    qq{SELECT queryId, desc,
              hasSpecific, hasSpecificStrong, hasSpecificCons,
              hasCofit, hasCofitCons
       FROM Query $whereClause},
    { Slice => {} });

print header,
    Utils::start_page("Fitness BLAST - $title"),
    '<div id="ntcontent">',
    h2($title);

if ($index > 0) {
    my $newIndex = $index - $maxToShow;
    $newIndex = 0 if $newIndex < 0;
    print p(a({-href => "batch_genes.cgi?jobId=$jobId&search=$searchSpec&hypo=$onlyHypo&link=$linkSpec&index=$newIndex"},
              "Previous page"));
}
my $maxIndex = $index+$maxToShow;
$maxIndex = scalar(@$qgenes) if $maxIndex > scalar(@$qgenes);
if ($index > 0 || scalar(@$qgenes) > $maxToShow) {
    print p(sprintf("Results %d to %d of %d", $index+1, $maxIndex, scalar(@$qgenes) ));
}

# And show the actual table
my @headings = map th($_), qw{Id Description Specific Cofit};
my @trows = ( $cgi->Tr({ -valign => 'middle', -align => 'center' }, @headings) );
foreach my $i ($index..($maxIndex-1)) {
    my $q = $qgenes->[$i];
    my $specString = ($q->{hasSpecificCons} ? "Cons." :
                      ($q->{hasSpecific} ? "Yes" : ""));
    my $specCol = $q->{hasSpecificStrong} ? "red" : "black";
    $specString = qq{<font color="$specCol">$specString</font>};
    my $specTitle = ();
    if ($q->{hasSpecific}) {
        $specTitle = $q->{hasSpecificStrong} ? "Has a strong and specific phenotype" : "Has a specific phenotype";
        $specTitle .= " and it is conserved" if $q->{hasSpecificCons};
    }

    my $cofitString = ($q->{hasCofitCons} ? "Cons." : ($q->{hasCofit} ? "Yes" : ""));
    $cofitString = qq{<font color="black">$cofitString</font>};
    my $cofitTitle = ($q->{hasCofitCons} ? "Conserved cofit" : ($q->{hasCofit} ? "Cofit (but it is not conserved)" : ""));

    push @trows, $cgi->Tr(td(a({-href => "batch_query.cgi?jobId=$jobId&queryId=$q->{queryId}"}, $q->{queryId})),
                          td(encode_entities($q->{desc})),
                          td(a({-title => $specTitle}, $specString)),
                          td(a({-title => $cofitTitle}, $cofitString)));
}
print table({cellpadding => 3, cellspacing => 0}, @trows);

if ($index + $maxToShow < scalar(@$qgenes)) {
    my $newIndex = $index + $maxToShow;
    print p(a({-href => "batch_genes.cgi?jobId=$jobId&search=$searchSpec&hypo=$onlyHypo&link=$linkSpec&index=$newIndex"},
              "Next page"));
}

print p(a({-href => "batch_overview.cgi?jobId=$jobId"}, "Overview of $jobname"));

Utils::endHtml($cgi);
