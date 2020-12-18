/* fitblast.js -- javascript code for fitness BLAST

   Fill an HTML element with a short summary of the fitness BLAST results:
   fitblast_load_short(element_id, server_root, sequence);

   Or fill an HTML element with a table of fitness BLAST results: 
   fitblast_load_table(element_id, server_root, sequence);

   See fitblast_example.html

   Copyright (c) 2015 University of California
   Authors: Morgan Price
 */

var mincoverage = 0.75; // hits below this coverage threshold are ignored
var mincloseid = 80; // %identity to be considered a close hit

var minabsstrong  = 2.0;

// a hit is considered likely to be useful if either
// cofitness is above threshold or (both absfit and absT are above thresholds)
var mincofit = 0.75; // minimum %cofitness considered useful
var minabsfit = 1.0; // minimum |fitness| considered useful
var minabsT = 4.0; // minimum |t| considered useful

function fitblast_usefulhit(row) {
    return(row.maxcofit >= mincofit
	       || (row.minfit < -minabsfit && row.minT < -minabsT)
	       || (row.maxfit > minabsfit && row.maxT > minabsT));
}

function fitblast_fitness_url(server_root, orgId, locusId) {
    return(server_root + 'cgi-bin/singleFit.cgi?orgId=' + orgId + '&locusId=' + locusId);
}

function fitblast_cofit_url(server_root, orgId, locusId) {
    return(server_root + 'cgi-bin/cofit.cgi?orgId=' + orgId + '&locusId=' + locusId);
}

function fitblast_locus_with_link(row, server_root) {
    var locus = row.name === "" ? (row.sysName === "" ? row.locusId : row.sysName) : row.name;
    return '<A HREF="' + fitblast_fitness_url(server_root,row.orgId,row.locusId) + '">' + locus + '</A>';
}

function fitblast_short(row,server_root) {
    var hascofit = row.maxcofit >= mincofit;
    var hassig = (row.minfit < -minabsfit && row.minT < -minabsT) || (row.maxfit > minabsfit && row.maxT > minabsT);
    var hasstrong = (row.minfit < -minabsstrong && row.minT < -minabsT) || (row.maxfit > minabsstrong && row.maxT > minabsT);
    var sigstring = hasstrong ? "strong phenotype" : (hassig? "has phenotype" : (row.minfit !== "" ? "has data" : "no data"));
    if (row.minfit !== "") {
	sigstring = '<A TITLE="fitness ' + (+row.minfit).toFixed(1) + " to " + (+row.maxfit).toFixed(1)
	    + '" HREF="' + fitblast_fitness_url(server_root,row.orgId,row.locusId) + '">'
	    + sigstring + "</A>";
    }
    if (hascofit) sigstring = sigstring + ', <A TITLE="top cofitness ' + (+row.maxcofit).toFixed(2) + '" HREF="'
		      + fitblast_cofit_url(server_root,row.orgId,row.locusId) + '">cofit</A>';
    return (+row.identity).toFixed(0) + '% id. to '
	+ fitblast_locus_with_link(row,server_root)
	+ ' from ' + row.organism + ': ' + sigstring;
}

// Show a short string with 1 or 2 lines with results. Results include links to fitness data.
function fitblast_load_short(id, server_root, sequence) {
   fitblast_loader(id, server_root, sequence, function(d, data) {
    // limit outbound links to 8000 a.a. of sequence
    // (instead of implementing POST)
    var detailURL = server_root + 'cgi-bin/mySeqSearch.cgi?query=' + sequence.substr(0,8000);
    var closeRow = null; 	// closest hit, if above threshold
    var pheRow = null;	// closest useful hit
    var i;
    for (i = 0; i < data.length; i++) {
      row = data[i];
      if (row.coverage < mincoverage) { continue; }
      var useful = fitblast_usefulhit(row);
      if (!closeRow && row.identity >= mincloseid) {
        closeRow = row;
        if (useful) { break; }
      } else if (useful) {
        pheRow = row;
        break;
      }
    }
    var out = [];
    if (closeRow) { out.push(fitblast_short(closeRow,server_root)); }
    if (pheRow) { out.push(fitblast_short(pheRow,server_root)); }
    if (!closeRow && !pheRow) { out.push("No hits with phenotypes"); }
    d.html( out.join("<BR>") + " (<A HREF='" + detailURL + "' TITLE='View more hits'>more</A>)" );
  });
}

function fitblast_load_table(id, server_root, sequence) {
   fitblast_loader(id, server_root, sequence, function(d,data) {
     data = data.filter(function(d) { return d.coverage >= mincoverage });
     if (data.length == 0) {
       d.html("No hits with high coverage");
       return;
     }
     //else
     var table = d.append('table');
     table.classed("fitblast",true);
     var thead = table.append('thead');
     var tbody = table.append('tbody');
     var columns = ['Identity', 'Organism', 'Locus', 'Description', 'Fitness', 'Cofit'];
     thead.append('tr').selectAll('th').data(columns)
       .enter().append("th").text(function(column) { return column; });
       var rows = tbody.selectAll("tr").data(data).enter().append("tr");
       rows.selectAll("td").data(function(row) {
         var showCofit = "";
         if (row.maxcofit >= mincofit) {
           showCofit = '<A HREF="' + fitblast_cofit_url(server_root, row.orgId, row.locusId) + '">'
           + (+row.maxcofit).toFixed(2) + '</A>';
         }
         var showFit = "No data";
         if (row.minfit !== "") {
           var colLo = row.minfit < -2 ? "#0000FF" : (row.minfit < -1 ? "#7777FF" : "#666666");
           var colHi = row.maxfit > 2 ? "#0000FF" : (row.maxfit > 1 ? "#AAAA00" : "#666666");
           showFit = '<span style="color:' + colLo + '">' + (+row.minfit).toFixed(1) + '</span>'
             + '<span style="color:#666666"> <small>to</small> </span>'
             + '<span style="color:' + colHi + '">' + (+row.maxfit).toFixed(1) + '</span>';
           if (fitblast_usefulhit(row)) { showFit = "<b>" + showFit + "</b>"; }
         }
         var show = { Organism: row.organism,
           Locus: fitblast_locus_with_link(row, server_root),
           Identity: row.identity + "%",
           Description: "<small>" + row.description + "</small>",
           Coverage: (+row.coverage).toFixed(2),
           Fitness: showFit,
           Cofit: '&nbsp; ' + showCofit
         };
         return columns.map(function(column) { return { column: column, value: show[column] } });
       }) //end computing columns for each row
       .enter().append("td").html(function(d) { return d.value; });
   });
}

// handler is given a document object to fill and the parsed hits
function fitblast_loader(id, server_root, sequence, handler) {
    var d = d3.select("#" + id);
    if (!d) {
	console.log("No element named "+id);
	return;
    }
    d.html("<small>loading...</small>");

    // for the main results, use POST to allow queries of over 8 kb or so
    // (otherwise apache considers it to long)
    var URL = server_root + "cgi-bin/seqservice.cgi";

    fetch(URL, {
      method: 'POST',
      headers: { 'Content-Type' : 'application/x-www-form-urlencoded' },
      body: "seq=" + sequence
    }).then(function(response) {
      if (!response.ok) {
        d.html("Cannot contact server");
      } else {
        return response.text();
      }
    }).then(function(text) {
      if (text) {
        data = d3.tsv.parse(text);
        if (!data || data.length < 1) {
          d.html("Empty response from server");
        } else if ("Error" in data[0]) {
          d.html(data[0].Error);
        } else {
          d.html("");
          handler(d, data);
        }
      }
   });
}
