#!/usr/bin/perl -w

#######################################################
## org.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page shows on overview  for a single 
# organism, with an overview of experiment groups nad
# links to specific phenotypes and downloads.
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $expGroup = $cgi->param('expGroup');
my $condition_1 = $cgi->param('condition_1'); 

my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};

# main table
# gather data and slice it into an array of hashes
my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup,
                                       COUNT(DISTINCT media || ":" || condition_1 || condition_2 || condition_3 || condition_4) AS nCond,
                                       COUNT(*) as nExp
                                       FROM Experiment WHERE orgId=?
                                       GROUP BY expGroup ORDER BY expGroup; },
                                    { Slice => {} },
                                    $orgId);


# write the title
my $title = scalar(@$cond) > 0 ? "$orginfo->{$orgId}{division}: $orginfo->{$orgId}{genome}" : "No experiments for this organism.";
my $start = Utils::start_page("$title");

print
    header,
    $start, '<div id="ntcontent">',
    h2($title);

#exit if no results
Utils::fail($cgi, "No experiments for this organism.") if @$cond == 0;

#create table
my @headings = qw{Group Conditions Experiments};
my @trows = ( Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
foreach my $row (@$cond) {
    push @trows, Tr({ -valign => 'top', -align => 'left' },
                    td([
                         $row->{expGroup} || "unrecorded", #group
                         $row->{nCond}, #conditions
                         a( { href => "exps.cgi?orgId=$orgId&expGroup=" . uri_escape($row->{expGroup}) },
                            $row->{nExp} ), #experiments
                       ]));
}

#display taxonomic information and link
if ((defined $orginfo->{$orgId}{taxonomyId}) && ($orginfo->{$orgId}{taxonomyId} ne "")) {
	print $cgi->p("NCBI taxonomy id:", $cgi->a({href => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$orginfo->{$orgId}{taxonomyId}"}, $orginfo->{$orgId}{taxonomyId}));
}

print
  h3("Fitness Experiments"),
  table({cellspacing => 0, style => "margin-left: 50px;", cellpadding => 3}, @trows),
  h3("Tools"),
  start_ul(),
  li(a( {href => "heatmap.cgi?orgId=$orgId"}, "Build heatmap")),
  li(a( {href => "keggmaplist.cgi?orgId=$orgId"}, "KEGG maps")),
  li(a( {href => "keggmap.cgi?mapId=01100&orgId=$orgId"}, "KEGG overview map")),
  li(a( {href => "seedsubsystemsOrg.cgi?orgId=$orgId"}, "SEED subsystems")),
  li(a( {href => "pathwaysOrg.cgi?orgId=$orgId"}, "MetaCyc pathways")),
  li(a( {href => "http://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=FitnessBrowser&gid=$orgId",
         title => "Find proteins in this genome that are relevant to a query term"},
        "Curated BLAST"), "for this genome"),
  li(a( {href => "http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?gdb=FitnessBrowser&gid=$orgId&set=aa",
         title => "Annotate amino acid biosynthesis pathways"},
        "GapMind for amino acid biosynthesis")),
  li(a( {href => "https://iseq.lbl.gov/genomes/genome/find?name=$orgId"},
        "ENIGMA genome browser")),
  end_ul();

autoflush STDOUT 1; # show preliminary results
print "\n";

# Overview of specific phenotypes
my $specN = $dbh->selectall_hashref(
    qq{ SELECT expGroup, count(DISTINCT condition_1) AS nCond, count(DISTINCT locusId) AS nGene
        FROM SpecificPhenotype JOIN Experiment USING (orgId,expName)
        WHERE orgId = ?
        GROUP BY expGroup },
    "expGroup", { Slice => {} }, $orgId);
if (scalar(keys %$specN) > 0) {
    my ($nSpecGene) = $dbh->selectrow_array(
        "SELECT COUNT(DISTINCT locusId) FROM SpecificPhenotype WHERE orgId=?",
        {}, $orgId);
    my ($nConsSpec) = $dbh->selectrow_array(
        "SELECT COUNT(DISTINCT locusId) FROM SpecOG WHERE orgId=? AND nInOG > 1",
        {}, $orgId);
    print
        h3("Specific Phenotypes"),
        p("$nSpecGene genes with specific phenotypes,  and $nConsSpec with conserved-specific phenotypes.");
    
    my @trows = ();
    push @trows, $cgi->Tr({-valign => "top"}, $cgi->th([ 'Group', '# Conditions', '# Genes' ]));
    foreach my $expGroup (sort keys %$specN) {
        my $row = $specN->{$expGroup};
        push @trows, $cgi->Tr({-valign => "top"}, $cgi->td([ 
                                  a({ href => "spec.cgi?orgId=$orgId&expGroup=".uri_escape($expGroup) }, $expGroup),
                                  $row->{nCond},
                                  $row->{nGene} ]));
    }
    print table({cellspacing => 0, style => "margin-left: 50px;", cellpadding => 3}, @trows);
}
print "\n";

my $reanno = $dbh->selectall_arrayref("SELECT * from Reannotation WHERE orgId=?",
                                      {}, $orgId);
if (scalar(@$reanno) > 0) {
    print
        h3("Updated Annotations"),
        a({-href => "orgReanno.cgi?orgId=$orgId"},
          "Updated annotations for " . scalar(@$reanno) . " genes");
}

my $sc = $dbh->selectall_arrayref("SELECT scaffoldId, length(sequence) AS len FROM ScaffoldSeq WHERE orgId=? ORDER BY len DESC",
                                  { Slice => {} }, $orgId);
my @scaffoldIds = map { $_->{scaffoldId} } @$sc;
my %scLabels = map { $_->{scaffoldId} => $_->{scaffoldId} . " (" . Utils::commify($_->{len}) . " nt)" } @$sc;

print h3("Scaffolds");
print
  start_form(-name => 'scaffolds', -method => 'GET', -action => 'genomeBrowse.cgi'),
  hidden('orgId'),
  p("Scaffold",
    popup_menu(-name => 'scaffoldId', -values => \@scaffoldIds, -labels => \%scLabels),
    "at",
    textfield( -name      => 'begin', -size => 8, -maxlength => 8 ),
    submit( -style => 'display: inline; text-align: left;', -name => 'Browse')),
  end_form;

# No COUNT DISTINCT in sqlite3 so use GROUP BY as workaround
my $loci = $dbh->selectcol_arrayref("SELECT locusId FROM ConservedCofit WHERE orgId = ? GROUP BY locusId",
                                          {}, $orgId);
my $nConsCofit = scalar(@$loci);
#file of experimental data - generated by createExpData.cgi
# my $data = `./createExpData.cgi orgId=$orgId`;
print h3("Downloads"),
  ul(li(a({href => "createFitData.cgi?orgId=$orgId"}, "Fitness values"), "(tab-delimited)"), 
     ul(li(a({href => "createFitData.cgi?orgId=${orgId}&t=1"}, "t scores"), "(tab-delimited)")),
     li(a({href => "createCofitData.cgi?orgId=$orgId"}, "Top cofitness for each gene"), "(tab-delimited)"),
     ul(li("includes $nConsCofit genes with conserved cofitness")),
     li(a({href => "spec.cgi?orgId=$orgId&download=1"}, "Specific phenotypes"), "(tab-delimited)"),
     li(a({href => "createExpData.cgi?orgId=$orgId"}, "Experiment meta-data"), "(tab-delimited)"),
     li(a({href => "downloadReanno.cgi?orgId=$orgId"}, "Reannotations"), "(tab-delimited; includes protein sequences)"),
     li(a({href => "orgSeqs.cgi?orgId=$orgId&type=nt"}, "Genome sequence"), "(fasta)"),
     li(a({href => "orgSeqs.cgi?orgId=$orgId"}, "Protein sequences"), "(fasta)"),
     li(a({href => "orgGenes.cgi?orgId=$orgId"}, "List of genes"), "(tab-delimited)")
    );
print "\n";

my $numGenes = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM Gene WHERE orgId = ?;}, undef, $orgId);
my $numData = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM GeneFitness WHERE orgId = ?;}, undef, $orgId);

print p("Fitness data for $numData of $numGenes genes.");

$dbh->disconnect();
Utils::endHtml($cgi);
