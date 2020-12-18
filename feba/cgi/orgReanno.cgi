#!/usr/bin/perl -w

#######################################################
## orgReanno.cgi
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################

# This page shows a list of reannotations for a single organism
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless exists $orginfo->{$orgId};

my $reanno = $dbh->selectall_arrayref(qq{ SELECT * from Reannotation JOIN Gene USING (orgId,locusId)
                                          WHERE orgId=?
                                          ORDER BY locusId },
                                      { Slice => {} }, $orgId);

my $title = "Updated Annotations for $orginfo->{$orgId}{genome}";
my $start = Utils::start_page("$title");

print
    header,
    $start, '<div id="ntcontent">',
    h2("Updated annotations for",
       a({-href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome}));

Utils::fail($cgi, "No genes in this organism have updated (or confirmed) annotations.") if @$reanno == 0;

print p(scalar(@$reanno), "genes with updated (or confirmed) annotations:");

foreach my $row (@$reanno) {
    my $showId = $row->{sysName} || $row->{locusId};
    my $locusId = $row->{locusId};
    $showId .= " $row->{gene}" if $row->{gene} ne "";
    my ($seed_desc) = $dbh->selectrow_array("SELECT seed_desc FROM SEEDAnnotation WHERE orgId = ? AND locusId = ?",
                                            {}, $orgId, $locusId);
    my $kegg_descs = $dbh->selectcol_arrayref("SELECT DISTINCT KgroupDesc.desc
                                               FROM BestHitKEGG JOIN KEGGMember USING (keggOrg,keggId)
                                               JOIN KgroupDesc USING (kgroup)
                                               WHERE orgId = ? AND locusId = ? AND desc <> ''",
                                              {}, $orgId, $locusId);
    my $kegg_desc = join("; ", @$kegg_descs);
    print p(a({-href => "singleFit.cgi?orgId=$orgId&locusId=$locusId"}, $showId),
            ":",
            $row->{new_annotation}, br(),
            small("Original description:", $row->{desc}), br(),
            small("SEED:", $seed_desc || "no annotation"), br(),
            small("KEGG:", $kegg_desc || "no annotation"), br(),
            small("Rationale:", $row->{comment}));
}

print p("Or download reannotations for",
        a({-href => "downloadReanno.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome}),
        "or for",
        a({-href => "downloadReanno.cgi"}, "all organisms"),
       "as tab-delimited tables");

$dbh->disconnect();
Utils::endHtml($cgi);
