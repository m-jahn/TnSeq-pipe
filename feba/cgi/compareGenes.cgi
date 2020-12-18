#!/usr/bin/perl -w
#######################################################
## compareGenes.cgi -- interactive scatterplot for two experiments
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Morgan Price
#######################################################
#
# Key parameeters: orgId, expName1 or query1, expName2 or query2
#	If query1 is set, expName1 is ignored, and similarly for query2
#	Cannot query on both simultaneously, however
# Optional: tsv -- use tsv=1 to fetch the data instead
#	outlier -- list outlying genes (xlow, xhigh, ylow, or yhigh)
#       with minabs -- minimum |abs| on selected axis.
# help -- 1 if on help/tutorial mode

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
print $cgi->header;

my $orgId = $cgi->param('orgId');
my $locus1 = $cgi->param('locus1');
my $locus2 = $cgi->param('locus2');
my $query1 = $cgi->param('query1');
my $query2 = $cgi->param('query2');
my $tsv = $cgi->param('tsv') ? 1 : 0;
my $help = $cgi->param('help') || "";
my $outlier = $cgi->param('outlier');
die "Must specify orgId" unless defined $orgId && $orgId ne "";
die "Must specify locus1 or query1"
  unless (defined $locus1 && $locus1 ne "")
  || (defined $query1 && $query1 ne "");
die "Must specify locus2 or query2"
  unless (defined $locus2 && $locus2 ne "")
  || (defined $query2 && $query2 ne "");
die "Cannot query both 1 and 2" if defined $query1 && defined $query2 && $query1 ne "" && $query2 ne "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless exists $orginfo->{$orgId};

my @geneCand = ();
my $choosing = undef;
if (defined $query1 && $query1 ne "") {
    my @genes = @{ Utils::matching_genes($dbh,$orgId,$query1) };
    Utils::fail($cgi, qq{No genes matching "$query1"}) if @genes == 0;
    @genes = grep { $_->{locusId} ne $locus2 } @genes if @genes > 1;
    if (@genes == 1) {
	$locus1 = $genes[0]{locusId};
    } else {
	@geneCand = @genes;
	$choosing = 1;
    }
} elsif (defined $query2 && $query2 ne "") {
    my @genes = @{ Utils::matching_genes($dbh,$orgId,$query2) };
    Utils::fail($cgi, qq{No genes matching "$query2"}) if @genes == 0;
    @genes = grep { $_->{locusId} ne $locus1 } @genes if @genes > 1;
    if (@genes == 1) {
	$locus2 = $genes[0]{locusId};
    } else {
	@geneCand = @genes;
	$choosing = 2;
    }
}

#choose if multiple
if (scalar(@geneCand) > 0) {
    die "Cannot use tsv mode with queries" if $tsv;
    # show table of these experiments
    my $locusConst = $choosing == 1 ? $locus2 : $locus1;
    my $geneConst = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
             {}, $orgId, $locusConst);
    die "Unknown locus: $locusConst" unless exists $geneConst->{locusId};
    my $notChoosing = $choosing == 1 ? 2 : 1;

    my @trows = ();
    my @headings = qw{&nbsp; LocusId Name SysName Description};
    push @trows, $cgi->Tr({-valign => 'top', -align => 'center'}, $cgi->th(\@headings));
    my $isFirst = 1;
    foreach my $gene (@geneCand) {
  my $checked = $isFirst ? "CHECKED" : "";
  push @trows, $cgi->Tr({-valign => 'top', -align => 'left'},
            $cgi->td([ qq{<input type="radio" name="locus$choosing" value="$gene->{locusId}" $checked >},
           $cgi->a({href => "geneOverview.cgi?orgId=$orgId&$gene->{locusId}"}, $gene->{locusId}),
           $gene->{gene}, $gene->{sysName}, $gene->{desc} ]));
  $isFirst = 0;
    }

  my $start = Utils::start_page("Select experiment to compare to");
    print $start, '<div id="ntcontent">',
  h2("Select gene in $orginfo->{$orgId}{genome}"),
  p("Selected gene will be compared to "
    . a( { href => "geneOverview.cgi?orgId=$orgId&locus=$locusConst" }, $geneConst->{gene} )
    . ": $geneConst->{desc}"),
  start_form(-name => 'input', -method => 'GET', -action => 'compareGenes.cgi'),
  hidden('orgId', $orgId),
  hidden("locus$notChoosing", $locusConst),
  table( {cellpadding => 3, cellspacing => 0}, @trows),
  submit('Go'),
  end_form;
    Utils::endHtml($cgi);
}


# orgId + locus > Gene Fitness > orgId + expName > Experiment

my $gene1 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
				   {}, $orgId, $locus1);
die "Unknown locus 1: $locus1" unless exists $gene1->{locusId};

my $gene2 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
				   {}, $orgId, $locus2);
die "Unknown locus 2: $locus2" unless exists $gene2->{locusId};

my $loc1 = $gene1->{locusId};
my $loc2 = $gene2->{locusId};
my $exps; # locusId => exps => attribute, with additional values x, y, tx, ty

my %fitVals;
my %tVals;
if ($tsv) {
    # fetch the data
    my $fit = $dbh->selectall_arrayref("SELECT * FROM GeneFitness WHERE orgId = ? AND locusId IN (?,?)",
               { Slice => {} }, $orgId, $loc1, $loc2);
    
    foreach my $row(@$fit) {
      $fitVals{$row->{locusId}}{$row->{expName}} = $row->{fit};
      $tVals{$row->{locusId}}{$row->{expName}} = $row->{t};
    }
    my $exps = $dbh->selectall_hashref("SELECT * FROM Experiment where orgId = ?", "expName", {}, $orgId);

    die "No experiments" unless scalar(keys %$exps) > 0;
    
    # my $found1 = 0;
    # my $found2 = 0;
    # foreach my $row (@$fit) {
    # 	my $locusId = $row->{locusId};
    # 	die "Unrecognized locus $locusId for org $orgId" unless exists $genes->{$locusId};
    # 	my $gene = $genes->{$locusId};
    # 	if ($row->{expName} eq $expName1) {
    # 	    $gene->{x} = $row->{fit};
    # 	    $gene->{tx} = $row->{t};
    # 	    $found1 = 1;
    # 	}
    # 	if ($row->{expName} eq $expName2) {
    # 	    $gene->{y} = $row->{fit};
    # 	    $gene->{ty} = $row->{t};
    # 	    $found2 = 1;
    # 	}
    # }
    # Utils::fail($cgi, "No fitness values for $expName1 in $orgId") unless $found1 > 0;
    # Utils::fail($cgi, "No fitness values for $expName2 in $orgId") unless $found2 > 0;

$dbh->disconnect();

    print join("\t", qw{expName expGroup expDesc x tx y ty})."\n";
    while (my ($expName, $expData) = each %$exps) {
    	next unless exists $fitVals{$locus1} && exists $fitVals{$locus2} && exists $exps->{$expName};
    	print join("\t", $expName, $exps->{$expName}->{expGroup}, $exps->{$expName}->{expDesc},
    		   $fitVals{$locus1}{$expName}, $tVals{$locus1}{$expName},$fitVals{$locus2}{$expName}, $tVals{$locus2}{$expName})."\n";
    }
    exit 0;

  }

# else interactive scatterplot

my $title = "Compare Genes for $orginfo->{$orgId}{genome}";
my $title2 = "Compare Genes for " . $cgi->a({href=>"org.cgi?orgId=$orgId"},$orginfo->{$orgId}{genome});
my $start = Utils::start_page("$title");

my $error = "";
if (Utils::gene_has_fitness($dbh,$orgId,$locus1) == 0) {
      $error = "Sorry! No fitness data for locus 1 (locus $locus1).";
    } if (Utils::gene_has_fitness($dbh,$orgId,$locus2) == 0) {
      $error = "Sorry! No fitness data for locus 2 (locus $locus2).";
    } 
    # else {
      # d3.select("#loading").html("Sorry! Cannot load data from " + tsvUrl + "<BR>Error: " + error);

my $showName1 = $gene1->{gene} || $gene1->{sysName} || $gene1->{locusId};
my $showName2 = $gene2->{gene} || $gene2->{sysName} || $gene2->{locusId};

my $showDesc1 = Utils::gene_link($dbh, $gene1, "desc", "domains.cgi");
my $showDesc2 = Utils::gene_link($dbh, $gene2, "desc", "domains.cgi");

print <<PART1
$start
<script src="../d3js/d3.min.js"></script>
<body style="padding-left: 1%">
<div id="ntcontent">

<H2>$title2</H2>

<P>
<i>x</i> axis <A HREF="geneOverview.cgi?orgId=$orgId&gene=$locus1">$showName1</A>: $showDesc1
<BR>
<i>y</i> axis <A HREF="geneOverview.cgi?orgId=$orgId&gene=$locus2">$showName2</A>: $showDesc2

  <p>$error</p>
PART1
;

if ($error ne "") {
  Utils::endHtml($cgi);
}


if ($help) {
    print qq[<div class="helpbox">
    <b><u>About this page:</u></b><BR><ul>
    <li>View the fitness values between two genes in this organism. </li>
    <li>To get to this page, search for any gene and click on the "Cofit" tab, then click on a cofitness value.</li> 
    <li>Change or flip the axes using the respective buttons and click on genes in the chart to add them to the table. Click on experiment names in the table to see more.</li>
    <li>To view the fitness heatmap for the two genes, click on the link below on the right hand side.</li>
    </ul></div>];
  }


print <<END
<div id="graphbox"> <div id="graphleft">


<div id="left"><!-- where SVG goes -->
<div id="loading"><!-- where status text goes -->
Please try another browser if this message remains
</div>
</div>
</div>

<div id="graphright">
<p>
<form method="get" action="compareGenes.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="locus2" value="$locus2" />
Change x axis: <input type="text" name="query1"  size="20" maxlength="100" />
<button type='submit'>Go</button>
</form>

<form method="get" action="compareGenes.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="locus1" value="$locus1" />
Change y axis: <input type="text" name="query2"  size="20" maxlength="100" />
<button type='submit'>Go</button>
</form>

<form method="get" action="compareGenes.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="locus1" value="$locus2" />
<input type="hidden" name="locus2" value="$locus1" />
<input type="submit" name="flip" value="Flip axes" />
</form>
</p>

  <P><b>Click on experiments to add them to the table:</b>

<TABLE id="genesel" cellspacing=0 cellpadding=3 >
<tr><th>Name</th><th>Group</th><th>Description</th><th>x</th><th>y</th><th>&nbsp;</th></tr>
</TABLE>
</P>

<P>
<form method="get" action="compareExps.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="locus1" value="$locus1" />
<input type="hidden" name="locus2" value="$locus2" />
<!--Or list outliers with
<select name="outlier">
   <option value="lowx">Low <i>x</i></option>
   <option value="lowy">Low <i>y</i></option>
   <option value="highx">High <i>x</i></option>
   <option value="highy">High <i>y</i></option>
</select>
and |fit| &gt; <select name="minabs" style="width: 60px;">
    <option value="1" selected>1.0 </option>
    <option value="1.5">1.5</option>
    <option value="2">2.0</option>
    <option value="2.5">2.5</option>
    <option value="3">3.0</option>
</select>
<input type="submit" name="submit" value="Go">-->
</form>

<P>
<A href="genesFit.cgi?orgId=$orgId&showAll=0&locusId=$locus1&locusId=$locus2">View heatmap for 2 genes</A>
  </P></P></div>
  <div style="clear: both;"></div></div>

<script>
var org = "$orgId";
var xLocus = "$gene1->{locusId}";
var yLocus = "$gene2->{locusId}";
var xName = "$showName1";
var yName = "$showName2";
var xDesc = "$gene1->{desc}";
var yDesc = "$gene2->{desc}";
var color = d3.scale.category10();

var margin = {top: 20, right: 130, bottom: 50, left: 75},
    width = 650 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width]);

var y = d3.scale.linear()
    .range([height, 0]);

//var color = d3.scale.category10();

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

var iSelected = 0; /* for color coding */
var selectColors = [ 'red', 'green', 'blue', 'magenta', 'brown', 'orange', 'darkturquoise' ];

var svg = d3.select("#left").append("svg")
    .attr("width",650)
    .attr("height",500)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

d3.select("#loading").html("Fetching data...");
var tsvUrl = "compareGenes.cgi?tsv=1&orgId=" + org + "&locus1=" + xLocus + "&locus2=" + yLocus;
d3.tsv(tsvUrl, function(error, data) {
  if (error || data.length == 0) {
      d3.select("#loading").html("Sorry! Cannot load data from " + tsvUrl + "<BR>Error: " + error);
      return;
  }
  d3.select("#loading").html("Formatting " + data.length + " genes...");
  data.forEach(function(d) {
    d.x = +d.x;
    d.y = +d.y;
    d.tx = +d.tx;
    d.ty = +d.ty;
  });

  var extentX = d3.extent(data, function(d) { return d.x; });
  var extentY = d3.extent(data, function(d) { return d.y; });
  var extentXY = d3.extent([ extentX[0], extentX[1], extentY[0], extentY[1] ]);
  x.domain(extentXY).nice();
  y.domain(extentXY).nice();

    svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("x", 350)
      .attr("y", 475)
      .style("text-anchor", "end")
      .text(xName);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("x", -80)
      .attr("y", -50)
      .style("text-anchor", "end")
      .text(yName);

  
  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(extentXY[0]))
       .attr("y2", y(extentXY[1]))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(0))
       .attr("y2", y(0))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.append("line")
       .attr("x1", x(0))
       .attr("x2", x(0))
       .attr("y1", y(extentXY[0]))
       .attr("y2", y(extentXY[1]))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

var cValue = function(d) { return d.expGroup;};
var tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0.0);

  svg.selectAll(".dot")
      .data(data)
    .enter().append("circle")
      .attr("class", "dot")
      .attr("r", 3)
      .attr("cx", function(d) { return x(d.x); })
      .attr("cy", function(d) { return y(d.y); })
      .style("fill", function(d) { return color(cValue(d)); })
      .style("z-index","1000")
      .on("click", dotClick)
      .on("mouseover", function(d) {
          tooltip.transition()
               .duration(200)
               .style("opacity", .9);
          tooltip.html(d.expDesc + ", " + d.expGroup + "<br/> (" + (+d.x).toFixed(1) 
          + ", " + (+d.y).toFixed(1)  + ")")
               .style("left", (d3.event.pageX + 5) + "px")
               .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseout", function(d) {
          tooltip.transition()
               .duration(500)
               .style("opacity", 0);
      });

  d3.select("#loading").html("");


  // legend = svg.append("g")
  //     .attr("class", "legend")
  //     .attr("transform","translate(50,30)")
  //     .style("font-size","12px")
  //     .call(d3.legend);

   // draw legend
  var legend = svg.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .attr("class", "legend")
      .style("font-size","10px")
      .attr("transform", function(d, i) { return "translate(100," + (i * 12 + 300) + ")"; });

  // draw legend colored rectangles
  legend.append("rect")
      .attr("x", width - 10)
      .attr("width", 10)
      .attr("height", 10)
      .style("fill", color);

  // draw legend text
  legend.append("text")
      .attr("x", width - 15)
      .attr("y", 5)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d;})

  /*
  var legend = svg.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .attr("class", "legend")
      .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });
  legend.append("rect").attr("x", width - 18).attr("width", 18).attr("height", 18).style("fill", color);
  legend.append("text").attr("x", width - 24).attr("y", 9).attr("dy", ".35em").style("text-anchor", "end").text(function(d) { return d; });
  */
});

function dotClick(d) {
    //var col = selectColors[(iSelected++) % selectColors.length];
    var col = color(d.expGroup);
    d3.select(this).attr("r",4).style("stroke","black");
    //.style("fill",col)
    columns = [ d.expName, d.expGroup, d.expDesc ];
    var tr = d3.select("#genesel").append("tr").attr("class","reset2").attr("valign","middle").style("color", col);
    //.style("color", col);
    var showId = d.expName;
    //d.sysName === "" ? d.locusId : d.sysName;
    var URL = "exp.cgi?orgId=" + org + "&expName=" + d.expName;
    var beginHref = "<A target='_blank' style='color: " + col + "' HREF='" + URL + "'>";
    tr.append("td").attr("class","expName").attr("expName",d.expName).html(beginHref + showId + "</A>");
    tr.append("td").html(d.expGroup);
    tr.append("td").html(d.expDesc);
    var hrefX = "strainTable.cgi?orgId=" + org + "&expName=" + d.expName + "&locusId="+xLocus;
    var hrefY = "strainTable.cgi?orgId=" + org + "&expName=" + d.expName + "&locusId="+yLocus;
    tr.append("td").html("<A style='color: " + col + "' TITLE='t = " + d.tx.toFixed(1) + "'"
                         + " HREF='"+hrefX+"'>" + d.x.toFixed(1) + "</A>");
    tr.append("td").html("<A style='color: " + col + "' TITLE='t = " + d.ty.toFixed(1) + "'"
                         + " HREF='"+hrefY+"'>" + d.y.toFixed(1) + "</A>");
    tr.append("td").html("<button type='button' onclick='removeRow(this)'>remove</button>");
}

function removeRow(a) {
    row = a.parentNode.parentNode; // href to td to row
    row.parentNode.removeChild(row);
}

function geneList() {
    var tds = document.getElementsByClassName("locusId");
    if (tds.length > 0) {
	var URL;
	if (tds.length == 1) {
	    URL = "myFitShow.cgi?orgId=" + org + "&gene=" + tds[0].getAttribute("locusId");
	} else {
            var i;
	    URL = "genesFit.cgi?orgId=" + org;
            for (i = 0; i < tds.length; i++ ) {
                URL += "&locusId=" + tds[i].getAttribute("locusId");
            }
	}
        window.open(URL);
   }
}

</script>

</div>
</body>
</html>
END
;

