#!/usr/bin/perl -w

#######################################################
## orgGroup.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# Shows the specific phenotypes for an organism and expGroup,
# or download all of the specific phenotypes in an organism
#
# Required CGI parameters:
# orgId -- which organism
#
# Optional CGI parameters, at least one of these must be specified:
# expGroup -- which group to show
# download -- set to 1 if in download mode

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use URI::Escape;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $dbh = Utils::get_dbh();
my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $expGroup = $cgi->param('expGroup');
my $download = $cgi->param('download');

if ($download) {
    die "Cannot specify expGroup with download" if $expGroup;
    my $spec = $dbh->selectall_arrayref(qq{SELECT expGroup, expName, condition_1, locusId, sysName, gene, desc, fit, t
                                           FROM SpecificPhenotype
                                           JOIN Experiment USING (orgId,expName)
                                           JOIN Gene USING (orgId,locusId)
                                           JOIN GeneFitness USING (orgId,locusId,expName)
                                           WHERE orgId = ?
                                           ORDER BY expGroup, condition_1, locusId, fit },
                                         { Slice => {} }, $orgId);
    my $cons = $dbh->selectall_arrayref(qq{SELECT locusId,expGroup,condition
                                             FROM SpecOG
                                             WHERE orgId = ? AND nInOG > 1 },
                                          { Slice => {} }, $orgId);
    my %cons = (); # expGroup => condition_1 => locusId => 1 if nInOG > 1
    foreach my $row (@$cons) {
        $cons{ $row->{expGroup} }{ $row->{condition} }{ $row->{locusId} } = 1;
    }
    print ("Content-Type:application/x-download\n");
    print "Content-Disposition: attachment; filename=specific_phenotypes_$orgId.txt\n\n";
    print join("\t", qw{expGroup expName condition_1 locusId sysName gene desc fit t conserved})."\n";
    foreach my $row (@$spec) {
        print join("\t", $row->{expGroup}, $row->{expName}, $row->{condition_1},
                   $row->{locusId}, $row->{sysName}, $row->{gene}, $row->{desc}, $row->{fit}, $row->{t},
                   exists $cons{ $row->{expGroup} }{ $row->{condition_1} }{ $row->{locusId} } ? "TRUE" : "FALSE"
            )."\n";
    }
    exit(0);
}
# else
my $style = Utils::get_style();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};
Utils::fail($cgi, "No experiment group specified") unless $expGroup;
my $title = "Specific Phenotypes for $orginfo->{$orgId}{genome} in $expGroup experiments";
my $title2 = "Specific Phenotypes for " 
    . a({href=>"org.cgi?orgId=$orgId"},$orginfo->{$orgId}{genome})
    . " in $expGroup experiments";
print
    header,
    Utils::start_page($title), '<div id="ntcontent">',
    h2($title2);

my $spec = $dbh->selectall_arrayref(
    qq{SELECT condition, locusId, sysName, gene, orgId, desc,
              minFit, maxFit, nInOG
       FROM SpecOG
       JOIN Gene USING (orgId,locusId)
       WHERE orgId = ? AND expGroup = ?
       ORDER BY condition, locusId },
    { Slice => {} }, $orgId, $expGroup);

if (@$spec == 0) {
    print "No specific phenotypes for $expGroup experiments";
}
my @headings = ("Condition", "Gene", "Description",
                a({ -title => "The strongest fitness among experiments with the specific phenotype" }, "Fitness"),
                a({ -title => "Does an ortholog have a specific phenotype in this condition?"}, "Cons.")
    );
my @trows = ();
push @trows, Tr({ -valign => 'top', -align => 'left' }, map { th($_) } \@headings);

my $lastCond = "";
foreach my $row (@$spec) {
    my $cond = $row->{condition};
    my $first = $cond ne $lastCond;
    $lastCond = $cond;
    my $cmpURL = "orthFit.cgi?orgId=$orgId&locusId=$row->{locusId}"
        . "&expGroup=" . uri_escape($expGroup)
        . "&condition1=" . uri_escape($cond);
    my $fitAttr = { -href => $cmpURL,
                    -style => "color:rgb(0,0,0)",
                    -title => $row->{minFit} != $row->{maxFit} ? sprintf("%.1f to %.1f", $row->{minFit}, $row->{maxFit}) : "no replicates" };
    my $allURL = "orthCond.cgi?expGroup=".uri_escape($expGroup)."&condition1=".uri_escape($cond);
    my $condParams = {href => $allURL, title=>"Compare specific phenotypes across organisms"};
    $condParams->{name} = $cond if $first;
    my $showFit = abs($row->{minFit}) > abs($row->{maxFit}) ? $row->{minFit} : $row->{maxFit};
    push @trows, Tr({ -valign => 'top', -align => 'left' },
                    td(a($condParams, $cond)),
                    td(Utils::gene_link($dbh, $row, "name", "myFitShow.cgi")),
                    td(Utils::gene_link($dbh, $row, "desc", "domains.cgi")),
                    td({ -bgcolor => Utils::fitcolor($showFit) },
                       a($fitAttr,
                         ( ($row->{minFit} < 0 && $row->{maxFit} < 0) || ($row->{minFit} > 0 && $row->{maxFit} > 0) ?
                           sprintf("%.1f", $showFit) : "variable" ))),
                    td($row->{nInOG} > 1 ? "Yes" : "&nbsp;"));
}
print table({cellspacing => 0, cellpadding => 3}, @trows); # style=>"margin-left: 0px",

$dbh->disconnect();
Utils::endHtml($cgi);

