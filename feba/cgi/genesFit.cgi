#!/usr/bin/perl -w
#######################################################
## genesFit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu), and 
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# locusId -- 1 or more locusIds
# Optional CGI parameters:
# addgene -- a gene to add (locusId, name, sysName)
# showAll -- 1 if showing all fitness values instead of just the most extreme ones
# around -- show neighboring genes (only with just 1 locusId, no addgene)
# help -- 1 if on help/tutorial mode

use strict;

use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util 'sum';
use HTML::Entities;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $orgId = $cgi->param('orgId') || die "No orgId found";
my @locusIds = $cgi->param('locusId');
if (@locusIds == 0) {
    print $cgi->header;
    Utils::fail($cgi,"no genes to show");
}
my $showAll = $cgi->param('showAll') ? 1 : 0;
my $around = $cgi->param('around') || 0;
my $addgene = $cgi->param('addgene');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";
my $expinfo = Utils::expinfo($dbh,$orgId);
my $help = $cgi->param('help') || "";

my $addgene_error = undef;
if ($addgene) {
    $addgene =~ s/[ \t]+$//;
    $addgene =~ s/^[ \t]+//;
    if ($addgene !~ m/^[A-Za-z90-9_.-]*$/) {
	$addgene_error = "Invalid gene to add";
    } else {
	my ($locusId) = $dbh->selectrow_array(qq{ SELECT locusId FROM Gene WHERE orgId = ?
                                                  AND (locusId = ? OR sysName = ? COLLATE NOCASE OR gene = ? COLLATE NOCASE) LIMIT 1 },
					      {}, $orgId, $addgene, $addgene, $addgene);
	if (defined $locusId) {
	    if (sum(map { $_ eq $locusId } @locusIds) > 0) {
		$addgene_error = qq{Gene "$addgene" (locus $locusId) is already included};
	    } else {
		push @locusIds, $locusId;
	    }
	} else {
	    $addgene_error = qq{Cannot find gene "$addgene" in $orginfo->{$orgId}{genome}};
	}
    }
}


my $centralId; # the focal gene
my $centralShow; # how to show its id
my %spacingDesc = (); # locusId => spacing description
my $type;
if ($around) {
    die "Cannot specify around with multiple locusIds or with addgene" unless @locusIds == 1;
    die "Invalid around = $around" unless $around =~ m/^\d+$/;
    $centralId = $locusIds[0];
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
				       {}, $orgId, $centralId);
    $centralShow = $gene->{sysName} || $gene->{locusId};
    $type = $gene->{type};
    my $scgenes = $dbh->selectall_arrayref("SELECT * from Gene where orgId = ? AND scaffoldId = ? ORDER BY begin",
					   { Slice => {} }, $orgId, $gene->{scaffoldId});
    die "Cannot find genes for $gene->{scaffoldId}" unless scalar(@$scgenes) > 0;
    my ($iCentral) = grep { $scgenes->[$_]{locusId} eq $centralId } (0..(scalar(@$scgenes)-1));
    die if !defined $iCentral;
    my $i1 = $iCentral - $around;
    $i1 = 0 if $i1 < 0;
    my $i2 = $iCentral + $around;
    $i2 = scalar(@$scgenes)-1 if $i2 >= scalar(@$scgenes);
    @locusIds = map { $scgenes->[$_]{locusId} } ($i1..$i2);
    foreach my $i ($i1..$i2) {
	my $leftsep = $i == 0 ? "" : $scgenes->[$i]{begin} - $scgenes->[$i-1]{end};
	my $rightsep = $i == scalar(@$scgenes)-1 ? "" : $scgenes->[$i+1]{begin} - $scgenes->[$i]{end};
	my $arrow = $scgenes->[$i]{strand} eq "+" ? "&#8594;" : "&#8592;"; # rightarrow or leftarrow
	$spacingDesc{ $scgenes->[$i]{locusId} } = join(" ",$leftsep,$arrow,$rightsep);
    };
}

my @genes = ();
foreach my $locusId (@locusIds) {
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
				       {}, $orgId, $locusId);
    die "No such locus $locusId in org $orgId" if !defined $gene->{locusId};
    # expName => "fit" => fitness value
    $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
                                               WHERE orgId = ? AND locusId = ?},
					   "expName", {}, $orgId, $locusId);
    foreach my $expName (keys %{ $gene->{fit} }) {
	die "No such experiment: $expName" unless exists $expinfo->{$expName};
    }
    $gene->{nExps} = scalar(keys %{ $gene->{fit} });
    push @genes, $gene;
}

my $pageTitle = $around ? "Fitness for $centralShow and surrounding genes in $genome"
    : "Fitness for " . scalar(@genes) . " genes in $genome";
my $start = Utils::start_page("$pageTitle");
my $tabs;
if ($around) {
    $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$centralId,0,$type,"nearby");
} else {
    $tabs = qq[<div id="ntcontent">];
}

print $cgi->header;
print $start, $tabs;

my @exps = values %$expinfo;
my $nTotValues = scalar(@exps);
if ($showAll) {
    @exps = sort Utils::CompareExperiments @exps;
} else {
    my $limitRows = 30;
    # for each experiment, compute average |fit| and average(fit)
    while (my ($expName,$exp) = each %$expinfo) {
	my @values = map { $_->{fit}{$expName}{fit} } @genes;
	@values = grep { defined $_ } @values;
	my @absvalues = map abs($_), @values;
	if (@values > 0) {
	    $exp->{absavg} = sum(@absvalues) / scalar(@values);
	    $exp->{avg} =  sum(@values)/scalar(@values);
	} else {
	    $exp->{absavg} = 0;
	    $exp->{avg} = 0;
	}
    }
    if (scalar(@exps) > $limitRows) {
	@exps = sort { $b->{absavg} <=> $a->{absavg} } @exps;
	@exps = @exps[0..($limitRows-1)];
    }
    @exps = sort { $a->{avg} <=> $b->{avg} } @exps;
}

print h2("Fitness for " . scalar(@genes) . " genes in " . $cgi->a({href => "org.cgi?orgId=". $orginfo->{$orgId}->{orgId}}, "$genome"),);

if (sum(map { $_->{nExps} > 0 ? 1 : 0} @genes) == 0) {
    print
        h3("No fitness data for these genes"),
        "<ul>";
    foreach my $gene (@genes) {
        print li(a({href => "geneOverview.cgi?orgId=$orgId&gene=$gene->{locusId}"},
                   $gene->{sysName} || $gene->{locusId})
                 . " : " . encode_entities($gene->{desc}) );
    }
    print "</ul>";
    $dbh->disconnect();
    Utils::endHtml($cgi);
    exit(0);
}


# corner box
print
    qq[<div style="position: relative;"><div class="floatbox">],
    
    start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
    "<P>Add gene: ",
    hidden( 'orgId', $orgId ),
    hidden( 'showAll', $showAll );
foreach my $locusId (@locusIds) { # avoid CGI sticky oddities
    print qq{<input type="hidden" name="locusId" value="$locusId" />\n};
}
print
    textfield( -name => 'addgene', -default => "", -override => 1, -size => 20, -maxLength => 100 ),
    # submit('Go'),
    "<button type='submit'>Go</button>",
    end_form;

print p,
    start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
    #"How many surrounding genes: ",
    "Number of genes on each side: ",
    hidden('orgId', $orgId),
    hidden('locusId', $centralId),
    hidden('showAll', $showAll),
    popup_menu(-name => 'around', -values => [ 1,2,3,4,5 ], -default => $around, -style => "min-width: 45pt"),
    "&nbsp;",
    # submit('Go'),
    "<button type='submit'>Go</button>",
    end_form
       if $around;
print "</P></div></div>";


if ($help) {
        print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>Show the experiments with the strongest phenotypes for the specified genes. (You can also choose to view all phenotypes with the link below.)</li>
        <li>To get to this page, search for any gene, visit the "Fitness" tab, and add genes to compare to.
        <li>Add more genes to compare to using the "Add gene" box at the top (try "uvrD").</li>
        <li>Or click on a "plot" link at the bottom to see a scatterplot.</li>
        <li>Hover over blue links for more information.</li>
        </ul></div>];
    }


print $cgi->h3($addgene_error) if defined $addgene_error;
if ($showAll) {
    print $cgi->p("All " . scalar(@exps) . " experiments, sorted by group and condition");
} else {
    print $cgi->p("Top " . scalar(@exps) . " experiments (either direction), sorted by average fitness");
}
print "Or view ";

my $locusSpec = $around ? "locusId=$centralId" : join("&", map {"locusId=$_"} @locusIds);
if ($showAll) {
    print $cgi->a( { href => "genesFit.cgi?orgId=$orgId&$locusSpec&around=$around" }, "strongest phenotypes" );
} else {
    print $cgi->a( { href => "genesFit.cgi?orgId=$orgId&$locusSpec&showAll=1&around=$around" },
                   "all $nTotValues experiments" );
}
print " or ",
  a({ -href => "heatmap.cgi?orgId=$orgId&" . join("&", map {"r=$_"} @locusIds) },
    "choose conditions");
print " or try the ",
  a({ -href => "cmpbrowser.cgi?anchorOrg=$orgId&" . join("&", map {"anchorLoci=$_"} @locusIds) },
    "comparative fitness browser") if $around;

print "<br><br>";


print Utils::geneArrows(\@genes, $centralId, undef, undef) if $around;

my @trows = ();
my @headings = qw{Group Condition};
my @headings2 = ("", "");
foreach my $gene (@genes) {
    push @headings, Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi");
    push @headings2, $gene->{gene};
}
push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings));

if (sum(map { $_ ne "" } @headings2) > 0) {
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings2));
}

# if ($around) {
#     my @headings3 = ("", "");
#     foreach my $gene (@genes) {
# 	my $spac = $spacingDesc{$gene->{locusId}};
#         push @headings3, qq{<div title="Spacing on either side of gene and strandedness"><small>$spac</small></div>};
#     }
#     push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->td(\@headings3));
# }

foreach my $exp (@exps) {
    my @values = ();
    my $expName = $exp->{expName};
    push @values, $cgi->td($exp->{expGroup});
    push @values, $cgi->td($cgi->a({ 
				     title => "$expName: $exp->{expDescLong}",
				     href => "exp.cgi?orgId=$orgId&expName=$expName" },
				   $exp->{expDesc}));
    foreach my $gene (@genes) {
	my $showId = $gene->{sysName} || $gene->{locusId};
	my $fit = $gene->{fit}{$expName}{fit};
	my $t = $gene->{fit}{$expName}{t};
    my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$gene->{locusId}&expName=$expName";
        $strainUrl .= "&help=1" if $help;
	if (defined $fit) {
	    my $fitShow = sprintf("%.1f",$fit);
	    $t = sprintf("%.1f", $t);
	    push @values, $cgi->td({ -bgcolor => Utils::fitcolor($fit) },
                                   a({ -href => $strainUrl,
                                       -title => "$showId: t = $t",
                                       -style => "color:rgb(0,0,0)" },
                                     $fitShow) );
	} else {
	    push @values, $cgi->td({ -bgcolor => Utils::fitcolor($fit) }, "&nbsp;");
	}
    }
    push @trows, $cgi->Tr({-align=>'left',-valign=>'top', -style => "font-size: 70%" }, @values);
}

my $compare = "";
my $firstGene = $centralId || $genes[0]->{locusId};
if (@genes > 0) {
    my @footer = ("","");
    foreach my $gene (@genes) {
	my @others = grep { $_->{locusId} ne $gene->{locusId} } @genes;
	my $url = "genesFit.cgi?orgId=$orgId&" . join("&", map { "locusId=$_->{locusId}" } @others);

    if ($around and ($centralId ne $gene->{locusId})) {
        $compare = "compareGenes.cgi?orgId=$orgId&locus1=$centralId&locus2=" . $gene->{locusId};
        $compare = "<BR>" . $cgi->a({href=>$compare}, "plot");
    } elsif (!$around and ($firstGene ne $gene->{locusId}) or ($firstGene ne $gene->{locusId})) {
        $compare = "compareGenes.cgi?orgId=$orgId&locus1=$firstGene&locus2=" . $gene->{locusId};
        $compare = "<BR>" . $cgi->a({href=>$compare}, "plot");
    }
    else {
        $compare = "";
    }
    push @footer, $cgi->a( { href => $url }, "remove")
			       . "<BR>" . ($gene->{sysName} || $gene->{locusId}) . $compare ;
    }
    push @trows, $cgi->Tr( { -align=>'center', -valign=>'top' }, $cgi->td(\@footer));
}

print table( { cellspacing => 0, cellpadding => 3 }, @trows);

print qq[<br><BR>];


$dbh->disconnect();
Utils::endHtml($cgi);
