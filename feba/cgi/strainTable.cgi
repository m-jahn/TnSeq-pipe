#!/usr/bin/perl -w
#######################################################
## strainTable.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# and either
#	scaffoldId, begin, end -- the range of interest (limited to 50 kb)
# or
#	locusId -- which gene to show (by default, range is exactly the width of the gene)
# or
#	scaffoldId, object (see below). [Object can also be specified together with begin/end]
# Optional CGI parameters:
# expName -- which experiment(s) to show. (Can be more than one, or none.)
# addexp -- additional experiments (i.e. a setname or a condition)
# zoom -- in or out
# pan -- left or right
# object -- an additional object to show, specified by a colon-delimited set of arguments such as
#	b:1000:e:2000:s:1:n:name:d:popup
#	begin, end, and name must be specified, and begin should be less than end
#	object is shown as unstranded if strand is absent
#	use 1 for + strand, not "+", because of CGI encoding issues
#	d is shown as popup text and is optional
#	(unlike genomeBrowse.cgi, just one object is supported)

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use StrainFitness;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";

my @expNames = $cgi->param('expName');
my $scaffoldId = $cgi->param('scaffoldId');
my $begin = $cgi->param('begin');
my $end = $cgi->param('end');
my $locusSpec = $cgi->param('locusId');
my $locusSpecShow;
my $tsv = $cgi->param('tsv') || 0;
my $expName = $cgi->param('expName') || "";
my $debug = $cgi->param('debug') || "";
my $help = $cgi->param('help') || "";
my $objspec = $cgi->param('object') || "";

if (defined $locusSpec && $locusSpec ne "") {
    my $sysName;
    ($scaffoldId,$begin,$end,$sysName) = $dbh->selectrow_array(
        "SELECT scaffoldId,begin,end,sysName FROM Gene WHERE orgId = ? AND locusId = ?",
        {}, $orgId, $locusSpec);
    die "Unknown locus $locusSpec in $orgId" if !defined $end;
    $locusSpecShow = $sysName || $locusSpec;
    # my $widen = int(1 + 0.2 * ($end-$begin+1));
    my $widen = 1000;
    $begin -= $widen;
    $end += $widen;
} elsif (defined $scaffoldId && defined $begin && defined $end) {
    die "Invalid scaffold parameter" if $scaffoldId eq "";
    die "Invalid begin parameter" unless $begin =~ m/^-?\d+$/;
    die "Invalid end parameter" unless $end =~ m/^\d+$/;
}

my %object = ();
if ($objspec) {
  die "object is used by scaffold is not set" unless $scaffoldId ne "";
  my %keynames = ("b" => "begin", "e" => "end", "s" => "strand", "n" => "name", "d" => "desc");
  my @parts = split /:/, $objspec;
  while (@parts > 0) {
    my $key = shift @parts;
    die "Unknown key $key in object argument\n" unless exists $keynames{$key};
    my $value = shift @parts;
    die "Wrong number of fields in $objspec\n" unless defined $value;
    $object{$keynames{$key}} = $value;
  }
  $object{object} = 1 ; # so that geneArrows() knows it is an object
  die "Invalid object $objspec\n"
    unless exists $object{begin} && exists $object{end} && exists $object{name};
  die "begin must be less than end" unless $object{begin} < $object{end};
  # To build the URL for getNtSeq.cgi, may need to reverse begin and end
  my ($beginL, $endL) = ($object{begin}, $object{end});
  ($beginL,$endL) = ($endL,$beginL)
    if exists $object{strand} && ($object{strand} eq "-" || $object{strand} eq "-1");
  # this link is currently always on the + strand, should fix
  $object{URL} = "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$beginL&end=$endL";

  if (!defined $begin || !defined $end) {
    $begin = $object{begin} - 500;
    $end = $object{begin} + 500;
    if ($begin - $end + 1 < 4000) {
      my $center = int(($begin+$end)/2);
      $begin = $center - 2000;
      $end = $center + 2000;
    }
  }
}

my $zoom = $cgi->param('zoom') || "";
my $initwidth = $end - $begin + 1;
if ($zoom eq "in") {
    $begin += 0.2 * $initwidth;
    $end -= 0.2 * $initwidth;
} elsif ($zoom eq "out") {
    $begin -= 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
my $pan = $cgi->param('pan') || "";
if ($pan eq "left") {
    
    $begin -= 0.4 * $initwidth;
    $end -= 0.4 * $initwidth;
} elsif ($pan eq "right") {
    $begin += 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
$begin = int($begin);
$end = int($end);

$end = $begin + 1 if $begin eq $end;
unless ($begin < $end) {
    print $cgi->header;
    Utils::fail($cgi, "Invalid begin/end $begin $end")
}
my $maxWidth = 40*1000;
if ($end - $begin >= $maxWidth) {
    print $cgi->header;
    Utils::fail($cgi, "Too wide, the limit is to show a range of " . Utils::commify($maxWidth) . " nucleotides");
}

my $addexp = $cgi->param('addexp');

# additional experiments?
if (defined $addexp && $addexp ne "") {
    my @expsNew = @{ Utils::matching_exps($dbh,$orgId,$addexp) };
    if (@expsNew == 0) {
        print header;
        Utils::fail($cgi, qq{No experiment matching "$addexp"});
    }
    # else
    my %oldExps = map { $_ => 1 } @expNames;
    push @expNames, grep {!exists $oldExps{$_} } map { $_->{expName} } @expsNew;
}

my $expinfo = Utils::expinfo($dbh,$orgId);
foreach my $expName (@expNames) {
    die "No such experiment: $expName" unless exists $expinfo->{$expName};
}

my $begComma = Utils::commify($begin);
my $endComma = Utils::commify($end);

#make tsv here? debate in printing first vs. running db commands
print header;

my $genes = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE orgId = ? AND scaffoldId = ? AND Gene.end >= ? AND Gene.begin <= ? ORDER by Gene.begin",
                               { Slice => {} }, $orgId, $scaffoldId, $begin, $end);
my %genesh = map {$_->{locusId} => $_} @$genes;

if ($tsv != 1) {
    my @expShow = ();
    foreach my $expName (@expNames) {
        push @expShow, a({href => "exp.cgi?orgId=$orgId&expName=$expName", title => $expName}, $expinfo->{$expName}{expDesc});
    }
    print
        Utils::start_page("Strain Fitness in $genome"),
        q{<div id="ntcontent">},
        h2("Strain Fitness in ",
           a({-href => "org.cgi?orgId=$orgId"}, "$genome"),
           defined $locusSpecShow ? "around " . a({-href => "singleFit.cgi?orgId=$orgId&locusId=$locusSpec"}, $locusSpecShow)
           : " at $scaffoldId: $begComma to $endComma"),
        p("Experiments: ",
          scalar(@expShow) > 0 ? join(", ", @expShow) : "none selected");

if ($help) {
        print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View the fitness of genes in various strains under a condition.</li>
        <li> Data on + strands are colored green, data on - strands are colored red, and data not affiliated with a gene are colored gray. Hover on points or blue links to see more information. Use the buttons to navigate the genome.</li>
        <li>To get to this page, click on the colored fitness boxes across the site (mostly on any pages relating to genes).</li> 
        </ul></div>];
    }

    print
        start_form(-name => 'input', -method => 'GET', -action => 'strainTable.cgi'),
        hidden( -name => 'orgId', -value => $orgId, -override => 1),
        hidden( -name => 'scaffoldId', -value => $scaffoldId, -override => 1),
        hidden( -name => 'begin', -value => $begin, -override => 1),
        hidden( -name => 'end', -value => $end, -override => 1),
        hidden( -name => 'object', -value => $objspec, -override => 1),
        join("\n", map { hidden( -name => 'expName', -value => $_, -override => 1) } @expNames),
        p({-class => "buttons", style=>"align:left; white-space:nowrap; line-height:40px;"}, "Add experiment(s): ",
           textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
           submit('Add','Add') ),
        p({-class => "buttons", style=>"max-width:500px; line-height:40px; white-space:nowrap;"},
          "Zoom:", submit('zoom','in'), submit('zoom','out'), "\tPan:", submit('pan','left'), submit('pan','right')),
        end_form;
    print 
        p(small(qq{Only strains with sufficient reads to estimate fitness are shown,
                   but the strain fitness values are still rather noisy and may be biased towards zero.
                   Strains near the edge of a gene are not shown as being associated with that
                   gene (they are in grey).}))
        if scalar(@expNames) > 0;
    

  if (@$genes == 0 && ! $objspec) {
      print "No genes in range.";
  } else {
      # sort @$genes;
      if ($debug) {
          foreach my $genea(@$genes) {
              print $genea->{begin} . "\t" . $genea->{end} . "\t | ";
          }
      }
      my @show = @$genes;
      push @show, \%object if keys(%object) > 0;
      @show = sort { $a->{begin} <=> $b->{begin} } @show;
      print Utils::geneArrows(\@show, $locusSpec, $begin, $end);
  }
}

my $tsvUrl = "strainTable.cgi?tsv=1&orgId=" . $orgId . "&scaffoldId=" . $scaffoldId . "&begin=" . $begin . "&end=" . $end . "&" . join("&", map {"expName=$_"} @expNames); #"&expName=" + expName;


# should I add zoom in/out and pan left/right buttons??
my $rows = StrainFitness::GetStrainFitness("../cgi_data", $dbh, $orgId, $scaffoldId, $begin, $end);

if (@$rows == 0) {
    print "No fitness data for strains in range " . Utils::commify($begin) . " to " . Utils::commify($end) . "\n";
}
my @trows = (); # the table
# header row
my @headings = qw{Position Strand Gene Sysname};
push @headings, a({-title => "Fractional position within gene"}, "Fraction");
my @base_headings = @headings;
foreach my $expName (@expNames) {
    push @headings, a({-href => "exp.cgi?orgId=$orgId&expName=$expName", -title => $expName},
                      $expinfo->{$expName}{expDesc});
}
push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, th(\@headings));

# leave out gene if not a used strain
foreach my $row (@$rows) {
    $row->{locusId} = "" unless $row->{used} eq "TRUE";
}
my %locusIds = map { $_->{locusId} => 1 } @$rows;
my %genes = (); # locusId => row
foreach my $locusId (keys %locusIds) {
    next if $locusId eq "";
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                       {}, $orgId, $locusId);
    die "Unknown locusId $locusId" unless exists $gene->{locusId};
    $genes{$locusId} = $gene;
}

# give up if no experiments specified
if (scalar(@expNames)==0) {
    exit 0 if $tsv==1; # empty output if tab-delimited mode
    Utils::endHtml($cgi);
}

# add row of links for removing items
my @remove_row = map { td("") } @base_headings; # 
my $baseURL = "strainTable.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end";
$baseURL .= "&object=$objspec" if $objspec;
foreach my $expName (@expNames) {
    my @otherExps = grep { $_ ne $expName } @expNames;
    my @otherExpSpec = map { "expName=$_" } @otherExps;
    push @remove_row, td( a({ -title => "remove $expName : $expinfo->{$expName}{expDesc}",
                              -href => join("&", $baseURL, @otherExpSpec) },
                            "remove") );
}
push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @remove_row);

my @avgFits = ();
foreach my $row (@$rows) {
    my $locusId = $row->{locusId};
    my $locusShow = "";
    my $gene = undef;
    if ($locusId ne "") {
        $gene = $genes{$locusId};
        $locusShow = $gene->{sysName} || $gene->{locusId};
    }
    my @row = ( a({-title => "barcode $row->{barcode}"}, Utils::commify($row->{pos})),
                $row->{strand},
                $genesh{$row->{locusId}}{gene},
                $locusId eq "" ? "" : a({-title => $gene->{desc}, -href => "singleFit.cgi?orgId=$orgId&locusId=$locusId"},
                                        $locusShow),
                $locusId eq "" ? "" : sprintf("%.2f",
                                              ($row->{pos} - $gene->{begin}) / ($gene->{end} - $gene->{begin} + 1))
        );
    @row = map { td($_) } @row;
    my $totalFit = 0; #gather the total for averaging
    my $ind = 0; #gather number of entries 
    foreach my $expName (@expNames) {
        my $fit = $row->{ $expName };
        $totalFit += $fit;
        $ind += 1;
        push @row, td( { -bgcolor => Utils::fitcolor($fit) }, sprintf("%.1f",$fit));
    }
    # above we give up if no experiments are specified, so this no longer causes divide by 0
    push @avgFits, $totalFit/$ind;

    push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @row);
}

if ($tsv == 1) { # tab delimited values, not a page

    print join("\t", qw{position strand gene locusId fit})."\n";
    my $ind = 0;
    foreach my $row (@$rows) {
        # next unless exists $gene->{x} && exists $gene->{y};
        my $displayName = $genesh{$row->{locusId}}{gene} || $genesh{$row->{locusId}}{sysName} || "";
        print join("\t", $row->{pos}, $row->{strand}, $displayName, $row->{locusId}, $avgFits[$ind])."\n";
        $ind += 1;
    }
    exit 0;
}

print <<END
<script src="../d3js/d3.min.js"></script>

<P>
<!--<i>x</i> axis: Position
<BR>
<i>y</i> axis: Average Strain Fitness-->

<TABLE width=100% style="border: none;">
<TR class="reset">
<TD valign="top" align="center" style="border: none;"><!-- left column -->

<div id="left"><!-- where SVG goes -->
<div id="loading"><!-- where status text goes -->
Please try another browser if this message remains
</div>
</div>
</TD>
</TR>


</TABLE>
</P>

<script>
var org = "$orgId";
var scaffoldId = "$scaffoldId";
var begin = "$begin";
var end = "$end";

var xName = "Position (kb)";
var yName = "Average Strain Fitness";


var margin = {top: 20, right: 20, bottom: 50, left: 50},
    width = 850 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width]);

var y = d3.scale.linear()
    .range([height, 0]);

var color = d3.scale.category10();
var cValue = function(d) { return d.strand;};

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

var iSelected = 0; /* for color coding */
var selectColors = [ 'red', 'green', 'blue', 'magenta', 'brown', 'orange', 'darkturquoise' ]; //red for -, green for +

// var svg = d3.select("#left").append("svg")
//     .attr("width",900)
//     .attr("height",500)
//   .append("g")
//     .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var svg = d3.select("#left")
   .append("div")
   .classed("svg-container", true) //container class to make it responsive
   .append("svg")
   //responsive SVG needs these 2 attributes and no width and height attr
   .attr("preserveAspectRatio", "xMinYMin meet")
   .attr("viewBox", "0 0 900 500")
   //class to make it responsive
   .classed("svg-content-responsive", true)
   .append("g")
   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

d3.select("#loading").html("Fetching data...");
var tsvUrl = "$tsvUrl"; //"strainTable.cgi?tsv=1&orgId=" + org + "&scaffoldId=" + scaffoldId + "&begin=" + begin + "&end=" + end + "&expName=" + expName;
d3.tsv(tsvUrl, function(error, data) {
  if (error || data.length == 0) {
      d3.select("#loading").html("Cannot load data from " + tsvUrl + "<BR>Error: " + error);
      return;
  }
  d3.select("#loading").html("Formatting " + data.length + " genes...");
  data.forEach(function(d) {
    d.position = (+d.position)/1000;
    d.fit = +d.fit;
  });

  var begin2 = begin/1000;
  var end2 = end/1000;

  var extentX = d3.extent([begin2, end2]); //data, function(d) { return d.position; });
  var extentY = d3.extent(data, function(d) { return d.fit; });
  var extentXY = d3.extent([ extentX[0], extentX[1], extentY[0], extentY[1] ]);
  x.domain(extentX);//.nice();
  y.domain(extentY).nice();

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("x", 350)
      .attr("y", 500-25)
      .style("text-anchor", "end")
      .text(xName);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("x", -80)
      .attr("y", -35)
      .style("text-anchor", "end")
      .text(yName);

  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(0))
       .attr("y2", y(0))
       .style("stroke","darkgrey")
       .style("stroke-width",1);


var tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0.0);


  svg.selectAll(".dot")
      .data(data)
    .enter().append("circle")
    .filter(function(d) { return d.locusId != "" })
      .attr("class", "dot")
      .attr("r", 5)
      .attr("cx", function(d) { return x(d.position); })
      .attr("cy", function(d) { return y(d.fit); })
      .style("fill", function(d) { 
        if (d.strand == '-'){return "red"} 
        else {return "green"}
        ; })
      .on("mouseover", function(d) {
          tooltip.transition()
               .duration(200)
               .style("opacity", .9);
          tooltip.html((d.gene||d.sysName||d.locusId) + ", at position " + (+d.position)+ " on " + d.strand + " strand, with fitness " + (+d.fit).toFixed(1))
               .style("left", (d3.event.pageX + 5) + "px")
               .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseout", function(d) {
          tooltip.transition()
               .duration(500)
               .style("opacity", 0);
      });


  svg.selectAll("dot")
        .data(data)
      .enter().append("circle")
      .filter(function(d) { return d.locusId == "" })
        .style("fill", "gray")
        .attr("r", 3.5)
        .attr("cx", function(d) { return x(d.position); })
        .attr("cy", function(d) { return y(d.fit); })
        .on("mouseover", function(d) {
          tooltip.transition()
               .duration(200)
               .style("opacity", .9);
          tooltip.html("position " + (+d.position) + " on " + d.strand + " strand, with fitness " + (+d.fit).toFixed(1))
               .style("left", (d3.event.pageX + 5) + "px")
               .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseout", function(d) {
          tooltip.transition()
               .duration(500)
               .style("opacity", 0);
      });

  d3.select("#loading").html("");

});

</script>

END
;
    
print h3("Per-strain Table"), small(table({ cellspacing => 0, cellpadding => 3, }, @trows));

print p("Or see this region's",
        a({ -href => "getNtSeq.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end" },
          "nucleotide sequence"));

$dbh->disconnect();
Utils::endHtml($cgi);
