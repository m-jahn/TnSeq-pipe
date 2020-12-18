#!/usr/bin/perl -w
#######################################################
## batch_query.cgi -- show best hits for a gene from a batch comparison
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# jobId -- the job identifier
# queryId -- the query protein
# Optional CGI parameter: mode=BLAST

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use HTML::Entities;
use IO::Handle; # for autoflush

use lib "../lib";
use Utils;
use Batch;

my $cgi=CGI->new;

my $jobId = $cgi->param('jobId');
my $queryId = $cgi->param('queryId');
my $mode = $cgi->param('mode');
die "Must specify jobId and queryId" unless $jobId ne "" && $queryId ne "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $bdb = Batch::get_batch_dbh($jobId);
my $jobname = Batch::get_job_name($jobId,$bdb);

my $qinfo = $bdb->selectrow_hashref("SELECT * FROM Query WHERE queryId=?", {}, $queryId);
die "Invalid queryId" unless $qinfo->{queryId} eq $queryId;

my $queryURL = undef;
$queryURL = "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$1" 
    if $queryId =~ m/^VIMSS(\d+)$/;

print header,
    Utils::start_page("Fitness BLAST - gene $queryId in $jobname"),
    '<div id="ntcontent">',
    h2($mode eq "BLAST" ? "All Hits" : "Best Hits",
       " for ",
       defined $queryURL ? a({-href => $queryURL}, $queryId) : $queryId,
       " from ", a({-href => "batch_overview.cgi?jobId=$jobId"}, $jobname)),
    p(b("Description: ", encode_entities($qinfo->{desc})));

if ($mode eq "BLAST") {
    my $aaseqSpec = "'" . $qinfo->{aaseq} . "'";
    print <<END
<script src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></script>
<script src="http://fit.genomics.lbl.gov/images/fitblast.js"></script>
<div id="fitblast"></div>
<script>
fitblast_load_table("fitblast", "http://fit.genomics.lbl.gov/", $aaseqSpec);
</script>
END
;
    print p("Or see", a({-href=>"batch_query.cgi?jobId=$jobId&queryId=$queryId"}, "best hits"));

} else { # regular mode not FastBLAST mode
    autoflush STDOUT 1; # so preliminary results appear
    my $besthits = $bdb->selectall_arrayref("SELECT * from BestHit WHERE queryId=? ORDER BY bits DESC",
                                            { Slice => {} }, $queryId);
    if (@$besthits == 0) {
        print p("Sorry, no good hits for this protein");
    } else {
        my @headings = (th("Orth"),th("Identity"),th("Organism"),th("Description"),
                        th({-colspan=>2},
                           a({-title=>"range of values"}, "Fitness")),
                        th("Specific"),
                        th(a({-title=>"Maximum cofitness"}, "Cofit")));
        my @trows = ( $cgi->Tr({ -valign => 'middle', -align => 'center' }, @headings) );
        foreach my $row (@$besthits) {
            my $orgId = $row->{orgId};
            my $locusId = $row->{locusId};
            my $orthSign = "?";
            $orthSign = "N" if ! $row->{isBBH};
            $orthSign = "Y" if $row->{coverage} >= 0.8 && $row->{identity} >= 30 && $row->{isBBH};
            my $orthText = sprintf("Coverage %.0f%% (minimum of both ways), %s",
                                   100 * $row->{coverage},
                                   $row->{bbh} ? "bidirectional best hit" : "not bidirectional");
            my $gene = $dbh->selectrow_hashref(qq"SELECT * from Gene where orgId=? AND locusId=?",
                                               {}, $orgId, $locusId);
            # Try to fail gracefully if $locus is invalid because of database changes
            my $locusShow = $gene->{sysName} || $locusId;
            my $locusDesc = $gene->{desc} || $locusId;
            my $fitLo = td("");
            my $fitHi = td("");
            my $maxCofit = "";
            my $spec_string = "";
            my $fitURL = "#";
            my $cofitURL = "#";
            if (defined $gene->{locusId}) {
                ($maxCofit) = $dbh->selectrow_array(
                    "SELECT cofit FROM Cofit where orgId = ? AND locusId = ? AND rank=1 LIMIT 1",
                    {}, $orgId, $locusId);
                $maxCofit = sprintf("%.2f", $maxCofit) if defined $maxCofit;
                my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(
                    qq{ SELECT min(fit), max(fit), min(t), max(t)
                        FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
                    {}, $orgId, $locusId);
                $fitURL = "singleFit.cgi?orgId=$orgId&locusId=$locusId";
                if (defined $minFit) {
                    $cofitURL = "cofit.cgi?orgId=$orgId&locusId=$locusId";
                    my $fitTitle = sprintf("t = %.1f to %.1f", $minT, $maxT);
                    $fitLo = td({-bgcolor => Utils::fitcolor($minFit), -style=>'text-align: center;' },
                               a( { -title => $fitTitle,
                                    -style => "color: rgb(0,0,0)",
                                    -href =>  $fitURL },
                                  sprintf("%.1f", $minFit)));
                    $fitHi = td({-bgcolor => Utils::fitcolor($maxFit), -style=>'text-align: center;' },
                               a( { -title => $fitTitle,
                                    -style => "color: rgb(0,0,0)",
                                    -href =>  $fitURL },
                                  sprintf("%.1f", $maxFit)));
                    my $spec = $dbh->selectall_arrayref(
                        "SELECT * from SpecOG WHERE orgId=? AND locusId=? ORDER BY minFit",
                        { Slice => {} }, $orgId, $locusId);
                    my @specs = ();
                    foreach my $row (@$spec) {
                        my $col = $row->{minFit} < 0 ? "darkblue" : "darkgoldenrod";
                        my $group = ($row->{expGroup} eq "carbon source" ? " (C)" :
                                     ($row->{expGroup} eq "nitrogen source" ? " (N)" : ""));
                        push @specs, span({-style => "color: $col", 
                                           -title => sprintf("fitness %.1f to %.1f", $row->{minFit}, $row->{maxFit}) },
                                           $row->{condition} . $group);
                    }
                    $spec_string = join(", ", @specs);
                }
            }
            push @trows, $cgi->Tr( td(a({-title => $orthText}, $orthSign)),
                                   td(sprintf("%.0f%%", $row->{identity})),
                                   td(small($orginfo->{$orgId}{genome})),
                                   td(small(a({-href => $fitURL,
                                               -title => $locusShow},
                                              encode_entities($locusDesc)))),
                                   $fitLo, $fitHi,
                                   td(small($spec_string)),
                                   td(a({-href=>$cofitURL}, $maxCofit)));
        }
        print table({cellpadding => 3, cellspacing => 0}, @trows);
    }
    print "\n"; # so preliminary results appear
    if ($qinfo->{hasCofitCons}) {
        # First, list the BBHs
        my %qOrth = (); # orgId => orthId
        my %qOrthBits = (); # orgId => bits
        my %qOrthId = (); # orgId => identity
        foreach my $row (@$besthits) {
            next unless $row->{isBBH};
            my $orgId = $row->{orgId};
            my $orthId = $row->{locusId};
            $qOrth{$orgId} = $orthId;
            $qOrthBits{$orgId} = $row->{bits};
            $qOrthId{$orgId} = $row->{identity};
        }

        # Then, fetch their conserved cofitness
        my %oCofit = (); # orgId => hitId (the hit of the BBH, not the BBH itself) => cofitness, for cons cofit pairs only
        my %bbhHits = (); # orgId => hitId => org2 => hit2
        while (my ($orgId,$orthId) = each %qOrth) {
            my $cons = $dbh->selectall_arrayref(qq{SELECT hitId,cofit,orth_orgId,orth_locusId,orth_hitId
                                                   FROM ConservedCofit WHERE orgId=? AND locusId=? },
                                                {}, $orgId, $orthId);
            foreach my $consrow (@$cons)  {
                my ($hitId,$cofit,$org2,$orth2,$hit2) = @$consrow;
                next unless $qOrth{$org2} eq $orth2; # indirect cases seem problematic
                $oCofit{$orgId}{$hitId} = $cofit;
                $bbhHits{$orgId}{$hitId}{$org2} = $hit2;
            }
        }

        if (scalar(keys %oCofit) > 0) {
            print
                h3("Conserved Cofitness"),
                p("Grouped by ", a({-href => "help.cgi#ortholog"}, "orthology")),
                "\n";

            my %shown = (); # orgId => hitId => 1 if this has been listed already
            # most-similar BBHs first
            my @orgs = sort { $qOrthBits{$b} <=> $qOrthBits{$a} } (keys %qOrthBits);
            my $nGroup = 0;
            my @rows = (); # group number or empty, orgId, BBH, hitId, cofitness
            foreach my $orgId (@orgs) {
                my $hitcofit = $oCofit{$orgId};
                my @hits = keys %$hitcofit;
                # Highest cofitness (in the closest organism) first
                @hits = sort {  $hitcofit->{$b} <=> $hitcofit->{$a} } @hits;
                foreach my $hitId (@hits) {
                    next if exists $shown{$orgId}{$hitId};
                    $shown{$orgId}{$hitId} = 1;
                    my @group = ();
                    $nGroup++;
                    push @rows, [ $nGroup, $orgId, $qOrth{$orgId}, $hitId, $oCofit{$orgId}{$hitId} ];
                    my $ohit = $bbhHits{$orgId}{$hitId};
                    my @orgsHit = sort { $qOrthBits{$b} <=> $qOrthBits{$a} } (keys %$ohit);
                    foreach my $org2 (@orgsHit) {
                        my $hit2 = $ohit->{$org2};
                        next if exists $shown{$org2}{$hit2};
                        $shown{$org2}{$hit2} = 1;
                        push @rows, [ "", $org2, $qOrth{$org2}, $hit2, $oCofit{$org2}{$hit2} ];
                    }
                }
            }
            my @trows = ();
            push @trows, $cgi->Tr(th("Group"),
                                  th("Organism"),
                                  th("Identity"),
                                  th({-colspan=>2}, "Ortholog"),
                                  th({-colspan=>2}, "Cofit with"),
                                  th("Cofit"));

            my $bgcolor = "";
            foreach my $row (@rows) {
                my ($group,$orgId,$orthId,$hitId,$cofit) = @$row;
                    $bgcolor = ($group % 2 ? "#FFFFFF" : "#DDDDDD")
                        if $group ne "";
                my $gene1 = $dbh->selectrow_hashref(qq"SELECT * from Gene where orgId=? AND locusId=?",
                                                    {}, $orgId, $orthId);
                my $gene2 = $dbh->selectrow_hashref(qq"SELECT * from Gene where orgId=? AND locusId=?",
                                                    {}, $orgId, $hitId);
                my $id1 = $gene1->{sysName} || $orthId;
                my $desc1 = $gene1->{desc} || $orthId;
                my $id2 = $gene2->{sysName} || $hitId;
                my $desc2 = $gene2->{desc} || $hitId;
                my $url1 = "singleFit.cgi?orgId=$orgId&locusId=$orthId";
                my $url2 = "singleFit.cgi?orgId=$orgId&locusId=$hitId";
                my $urlCofit = "cofitCons.cgi?orgId=$orgId&locusId=$orthId&hitId=$hitId";
                push @trows, $cgi->Tr({ -bgcolor => $bgcolor},
                                      td($group),
                                      td(small($orginfo->{$orgId}{genome})),
                                      td(sprintf("%.0f%%", $qOrthId{$orgId})),
                                      td(small(a({-href => $url1}, $id1))), 
                                      td(small(encode_entities($desc1))),
                                      td(small(a({-href => $url2}, $id2))),
                                      td(small(encode_entities($desc2))),
                                      td(a({-href => $urlCofit, -title => "compare cofitness"},
                                           sprintf("%.2f", $cofit))));
            }
            print table({cellpadding => 3, cellspacing => 0}, @trows), p("&nbsp;");
        }
    }

    print p("Or see",
            a({-href=>"batch_query.cgi?jobId=$jobId&queryId=$queryId&mode=BLAST"},
              "all homologs for this protein"));
}

Utils::endHtml($cgi);

