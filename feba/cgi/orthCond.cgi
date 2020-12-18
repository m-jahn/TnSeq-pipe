#!/usr/bin/perl -w
#######################################################
## orthCond.cgi -- compare specific phenotypes for a given Group & Condition_1 across all bugs
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Given a group and condition, all specific-sick genes,
# grouped into ad hoc ortholog groups
#
# Required CGI parameters: expGroup and condition1 (condition1 can be empty, but requires exact match)
# Optional: help -- 1 if on help/tutorial mode

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;

sub RowForGene($$$$$); # gene object, shade or not, which row (0 indexed) => the row
sub summaryRow($$$); # ortholog group, shade, first to include (1 indexed)
sub CompareOG($$); # a comparator for sorting OGs (lists of genes)
sub OrderGenesInOG; # a comparator for sorting genes within an OG

my $cgi=CGI->new;
my $style = Utils::get_style();

my $expGroup = $cgi->param('expGroup') || die "no expGroup parameter";
my $condition1 = $cgi->param('condition1');
die "no condition1 parameter" unless defined $condition1;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $debug = $cgi->param('debug');
my $start_time = gettimeofday() if $debug;
my $help = $cgi->param('help') || "";

my $speclist = $dbh->selectall_arrayref(qq{SELECT ogId,orgId,locusId,minFit,maxFit,minT,maxT,sysName,gene,desc
                                           FROM SpecOG JOIN Gene USING (orgId,locusId)
					   WHERE expGroup = ? AND condition = ?},
					   { Slice => {} }, $expGroup, $condition1);

my $title = "Specific Genes for $expGroup Experiments in $condition1 Across Organisms";
$title = "No Specific Genes" if (@$speclist == 0);
my $start = Utils::start_page("$title");

my $js =  q[<script type="text/javascript" src="../images/jquery-1.11.3.min.js"></script>
<script type="text/javascript">$(document).ready(function(){ 
    $('tr.header2').nextUntil('tr.header').hide(); // after each tr.header2, hide rows until the next header
    $('tr.header2').click(function(){
	// click on header2 to hide it and turn on elements until next header
        $(this).toggle();
        $(this).nextUntil('tr.header').css('display', function(i,v){
            return this.style.display === 'table-row' ? 'none' : 'table-row';
        });
    });
    $('tr.header3').click(function(){
        // click on header3 to hide it and turn off elements until next header. show the previous row (header2)
        $(this).toggle();
        $(this).closest('tr').prev().show();
        // $('tr.header2').show();
        $(this).nextUntil('tr.header').css('display', function(i,v){
            return this.style.display === 'table-row' ? 'none' : 'table-row';
        });
    });
    
});
    </script>];

print
    header, $start, $js, '<div id="ntcontent">',
    h2($title);

Utils::fail($cgi,"No genes with specific phenotypes were found for this condition") if (@$speclist == 0);

# Group the results by ogId and sort the OGs by how many members they have
my %og = (); # ogId => list of rows
foreach my $spec (@$speclist) {
    push @{ $og{$spec->{ogId}} }, $spec;
}
my @ogInOrder = sort { CompareOG($og{$a}, $og{$b}) } (keys %og);

my @headings = ['&nbsp', 'Organism', 'Gene', 'Name', 'Description', 'Fitness (Lower)', 'Fitness (Upper)'];  #Experiment
my @trows = ( $cgi->Tr({ -valign => 'middle', -align => 'center' }, map { th($_) } @headings ) );

my $shade = 0;
my $group = 1;
my $singletonHead = 0;

foreach my $ogId (@ogInOrder) {
    my @sorted = sort OrderGenesInOG @{ $og{$ogId} };
    my $row = 0;
    foreach my $gene (@sorted) {
        if (scalar(@sorted) == 1 and $singletonHead == 0) {
            push @trows, qq[<th colspan="8"><center>Singletons</center></th>];
            $singletonHead = 1;
        }
        if ($row == 3 and scalar(@sorted) > 4) {
            push @trows, summaryRow(\@sorted, $shade, 4);
            push @trows, qq[<tr class="header3"><th colspan="8"><center><span>Collapse -</span></center></th></tr>];
        }
        push @trows, RowForGene($dbh, $gene, $shade, $singletonHead == 1 ? "" : $group, $row);
        $row++;
    }
    $shade++;
    $group++;
}

if ($help) {
        print qq[<BR><BR><div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View genes with <A HREF="help.cgi#specific">specific and strong phenotypes</A> in $condition1 $expGroup, across all organisms.</li>
        <li>These genes are grouped by <A HREF="ortholog">orthology</A>, and larger ortholog groups are shown first.</li>
        <li>To get to this page, search for any experiment and click on the link at the bottom to view specific phenotypes across organisms.</li> 
        <li>Click on "+" to expand each section.</li>
        <li>Click on a fitness value to see all of the fitness values for that gene and its orthologs in $condition1 $expGroup.
        </ul></div>];
    }


print
    p("Genes with",
      a({ -href => "help.cgi#specific" },  "specific phenotypes"),
      "in $expGroup $condition1 are",
      font({ style => "color: rgb(120,120,120); font-weight: bold;"}, "grouped"),
      "by",
      a({ -href => "help.cgi#ortholog"}, "orthology")),
    table({cellpadding => 3, cellspacing => 0}, @trows);

    print "<BR>";

$dbh->disconnect();
Utils::endHtml($cgi);

sub summaryRow($$$) {
    my ($og, $shade, $firstToInclude) = @_; # firstToInclude is 1-based
    my $i = 1;
    my @row;
    my @showOrgs = ();

    foreach my $gene(@$og) {
        my $orgId = $gene->{orgId};
        my $genome = $orginfo->{$orgId}{genome};
        my $genomeShort = $genome;
        $genomeShort =~ s/^(.)\S+/$1./;
        my $locusId = $gene->{locusId};

        push @showOrgs, $cgi->a({href=>"orthFit.cgi?orgId=$orgId&locusId=$locusId"
                                 . "&expGroup=" . uri_escape($expGroup)
                                 . "&condition1=" . uri_escape($condition1),
                                 title=>"$gene->{desc}"},
                                $genomeShort)
            if $i >= $firstToInclude;
        $i++;
    }
    
    my $count = scalar(@$og) - ($firstToInclude-1);
    push @row, $cgi->Tr(
        {-class=>'header2', -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"},
        td(span({class=>"deco"}, a({title=>"Expand"}, '+'))),
        td({colspan=>"6"}, "$count more from " . join(", ", @showOrgs))
    );
    return @row;
}

# $row is 0-indexed; $group is the ortholog group number (1-indexed)
sub RowForGene($$$$$) {
    my ($dbh,$gene,$shade,$group,$row) = @_;
   
    # my $firstForGene = 1;
    my $orgId = $gene->{orgId};
    my $genome = $orginfo->{$orgId}{genome};
    my $genomeShort = $genome;
    $genomeShort =~ s/^(.)\S+/$1./;
    my $locusId = $gene->{locusId};

    my $rowLabel = "";
    my $collapse = "";
    if ($row == 0) {
        $rowLabel = a({style=>"color: #011B47", title=>"Ortholog group $group"},$group);
        $collapse = 'header';    }
    my $orthFitURI = "orthFit.cgi?orgId=$orgId&locusId=$locusId"
        . "&expGroup=" . uri_escape($expGroup)
        . "&condition1=" . uri_escape($condition1);
    $orthFitURI .= "&help=1" if $help;

    return $cgi->Tr( { -class=> $collapse, -valign => 'middle', -align => 'left',
                       -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF" },
                     td($rowLabel),
                     td( a({href => "spec.cgi?orgId=$gene->{orgId}&expGroup="
                            . uri_escape($expGroup) . "#" . $condition1 },
                           $genomeShort) ),
                     td( Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi") ),
                     td( $gene->{gene} ),
                     td( Utils::gene_link($dbh, $gene, "desc", "domains.cgi") ),
                     td( { -bgcolor => Utils::fitcolor($gene->{minFit}), -style=>'text-align: center;' },
                         a( { -title => sprintf("Click to compare (t = %.1f to %.1f)",$gene->{minT},$gene->{maxT}),
                              -style => "color: rgb(0,0,0)",
                              -onMouseOver=>"this.style.color='#CC0024'",
                              -onMouseOut=>"this.style.color='#000000'",
                              -href => $orthFitURI },
                            sprintf("%.1f",$gene->{minFit}) ) ),
                     td( { -bgcolor => Utils::fitcolor($gene->{maxFit}), -style=>'text-align: center;' },
                         a( { -title => sprintf("Click to compare (t = %.1f to %.1f)",$gene->{minT},$gene->{maxT}),
                              -style => "color: rgb(0,0,0)",
                              -onMouseOver=>"this.style.color='#CC0024'",
                              -onMouseOut=>"this.style.color='#000000'",
                              -href => $orthFitURI },
                            sprintf("%.1f",$gene->{maxFit}) ) )
                   );
}

# sort them so the ones with gene names come up first;
# else put E. coli first; else sort by genome name
sub OrderGenesInOG {
    my $cmp = ($b->{gene} ne "") <=> ($a->{gene} ne "");
    return $cmp if $cmp;
    $cmp = ($b->{orgId} eq "Keio") <=> ($a->{orgId} eq "Keio");
    return $cmp if $cmp;
    my $genomeA = $orginfo->{ $a->{orgId} }{genome};
    my $genomeB = $orginfo->{ $b->{orgId} }{genome};
    return $genomeA cmp $genomeB;
}

# larger OGs first. Break ties by the lowest genome name
sub CompareOG($$) {
    my ($la,$lb) = @_; # each is a list
    my $cmp = scalar(@$lb) <=> scalar(@$la);
    return $cmp if $cmp;
    my @genomesA = sort map { $orginfo->{ $_->{orgId} }{genome} } @$la;
    my @genomesB = sort map { $orginfo->{ $_->{orgId} }{genome} } @$lb;
    $cmp = $genomesA[0] cmp $genomesB[0];
    return $cmp if $cmp;
    $cmp = ($la->[0]{sysName} || $la->[0]{locusId}) cmp 
           ($lb->[0]{sysName} || $lb->[0]{locusId});
    return $cmp;
}
