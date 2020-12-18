#!/usr/bin/perl -w
#######################################################
## exp.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price, Victoria Lo
#######################################################
#
# Required parameters: orgId and expName, for organism and which experiment
# Optional parameter show:
# specific -- show specific phenotypes
# important or detrimental -- show top 200 genes either way
# quality -- show quality metrics
# (by default, shows all of these optoins)
# help -- 1 if on help/tutorial mode

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
sub CompoundToHTML($$$$);
sub MediaComponentHTML($$);

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
die "Invalid species name!" unless $orgId =~ m/^[A-Za-z0-9_-]*$/;
my $expName = $cgi->param('expName') || die "no experiment identifier\n";
$expName =~ s/^[ \t]+//;
$expName =~ s/^[ \t\r\n]+$//;
die "Invalid experiment name!" unless $expName =~ m/^[A-Za-z0-9_]*$/;
my $show = $cgi->param('show') || "";
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "No such orgId: $orgId" unless exists $orginfo->{$orgId};

my $exp = $dbh->selectrow_hashref("SELECT * FROM Experiment WHERE orgId = ? AND expName = ?", {}, $orgId, $expName);
die "No such experiment: $expName" unless defined $exp->{expName};

# Fetch specific phenotypes
my $spec = $dbh->selectall_arrayref(qq{SELECT * from SpecificPhenotype JOIN GeneFitness USING (orgId,locusId,expName)
                                       WHERE orgId = ? AND expName = ?},
				    { Slice => {} },
				    $orgId, $expName);
my $start = Utils::start_page("Experiment $expName for $orginfo->{$orgId}{genome}");
my $tabs = "";
if ($show eq "") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"overview");
} elsif ($show eq "specific") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"specific");
} elsif ($show eq "important"){
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"-gene");
} elsif ($show eq "detrimental") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"+gene");
} elsif ($show eq "quality") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"metrics");
}

print $start, $tabs,
    h2("Experiment $expName for ". $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
	# corner compare box
	qq[<div style="position: relative;"><div class="floatbox">],
    start_form(-name => 'input', -method => 'GET', -action => 'compareExps.cgi'),
    hidden('orgId', $orgId),
    hidden('expName2', $expName),
    "Compare to: ",
    textfield(-name => 'query1', -value => '', -size => 20, -maxlength => 100),
    # submit('Go'),
    "<button type='submit'>Go</button>",
    end_form,
    qq[</P></div></div>],
    h3($exp->{expDescLong});

if ($help) {
  print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View the genes that had strong and specific phenotypes in this experiment. </li>
        <li>To get to this page, search for any experiment and click on the "Specific" tab.</li> 
        <li>To compare to another experiment via scatterplot, add another experiment using the box above. (Try "cisplatin".)</li>
        <li>To make a comparative heatmap, check the genes of interest and click the "Heatmap" link at the bottom.</li>
        <li>For more about how we define a specific phenotype, see the <A HREF="help.cgi#specific">help page</A>.
        </ul></div>];
}

autoflush STDOUT 1; # so preliminary results appear
print "\n";

my @fit = (); # sorted list of fitness values to show
my $header = undef;
if ($show eq "") {
    my $mediaList = $dbh->selectall_arrayref("SELECT DISTINCT media FROM MediaComponents");
    my %mediaOrMix = map { $_->[0] => 1 } @$mediaList;
    my @cond = ();
    foreach my $i (1..4) {
      push @cond, CompoundToHTML( $exp->{"condition_".$i}, $exp->{"concentration_".$i}, $exp->{"units_".$i},
                                \%mediaOrMix );
    }
    @cond = grep { $_ ne "" } @cond;
    my $media = $exp->{media};
    $media .= " ($exp->{mediaStrength}x)" if $exp->{mediaStrength} != 1;
    $media .= " + " . join(" + ", @cond) if @cond > 0;
    $media .= ", pH=$exp->{pH}" if $exp->{pH} ne "";
    my @culture = ("Culturing: ". $exp->{mutantLibrary});
    push @culture, $exp->{vessel} if $exp->{vessel} ne "";
    push @culture, $exp->{aerobic} if $exp->{aerobic} ne "";
    push @culture, "at $exp->{temperature} (C)" if $exp->{temperature} ne "";
    push @culture, "shaken=$exp->{shaking}" if $exp->{shaking} ne "";
    push @culture, "($exp->{liquid})" if $exp->{liquid} ne "" && lc($exp->{liquid}) ne "liquid";

    my @pieces = ("Media: $media", join(", ",@culture));
    push @pieces, sprintf("Growth: about %.1f ", $exp->{nGenerations})
      . a({-href => "help.cgi#growth"}, "generations")
        if $exp->{nGenerations};
    push @pieces, "By: $exp->{person} on $exp->{dateStarted}";
    if ($exp->{pubId}) {
      my $pub = $dbh->selectrow_hashref("SELECT * from Publication WHERE pubId = ?",
                                        {}, $exp->{pubId});
      if (defined $pub) { # should always exist
        push @pieces, "Reference: " . a({ -href => $pub->{URL}, -title => $pub->{title} },
                                     $pub->{pubId});
      }
    }

    my $mediaComponents = $dbh->selectall_arrayref(qq{SELECT * from MediaComponents LEFT JOIN Compounds USING (compound)
                                                      WHERE media = ? },
                                                   { Slice => {} },
                                                   $exp->{media});
    if (@$mediaComponents > 0) {
      my $html = MediaComponentHTML($mediaComponents, $exp->{mediaStrength});
      $html .= " " . small("(final concentrations)") if $exp->{mediaStrength} != 1;
      push @pieces, "Media components: $html";
    }
    foreach my $i (1..4) {
      if (exists $mediaOrMix{$exp->{"condition_".$i}}
          && lc( $exp->{"units_".$i} ) eq "x"
          && $exp->{"concentration_".$i} =~ m/^\d+[.]?\d*$/) {
        my $comp = $dbh->selectall_arrayref(qq{SELECT * FROM MediaComponents LEFT JOIN Compounds USING (compound)
                                              WHERE media = ? },
                                            { Slice => {} },
                                            $exp->{condition_1});
        my $html = MediaComponentHTML($comp, $exp->{concentration_1});
        $html .= " " . small("(final concentrations)") if $exp->{"concentration_".$i} != 1;
        push @pieces, $exp->{"condition_".$i} . " " . $exp->{"concentration_".$i} . "x includes: $html";
      }
    }
    if ($exp->{growthPlate} ne "" && $exp->{growthWells} ne "") {
      push @pieces, "Growth plate: $exp->{growthPlate} $exp->{growthWells}";
    }
    print join("<BR>\n", @pieces)."\n";
} elsif ($show eq "specific") {
    $header = "Genes with " . a({href => "help.cgi#specific"}, "specific") . " phenotypes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * FROM SpecificPhenotype JOIN GeneFitness USING (orgId,expName,locusId)
                                          JOIN Gene USING (orgId,locusId)
					  WHERE orgId=? AND expName=? ORDER BY fit},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "important") {
    $header = "200 most important genes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * from GeneFitness JOIN GENE USING (orgId,locusId)
                                          WHERE orgId=? AND expName=?
                                          ORDER BY fit LIMIT 200},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "detrimental") {
    $header = "200 most detrimental genes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * from GeneFitness JOIN GENE USING (orgId,locusId)
                                          WHERE orgId=? AND expName=?
                                          ORDER BY fit DESC LIMIT 200},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "quality") {
    $header = "Quality Metrics:";
}

my %cons = (); # locusId => nInOG if it has a specific phenotype in this condition
if ($exp->{condition_1} ne "") {
    my $specOG = $dbh->selectall_arrayref("SELECT locusId,nInOG from SpecOG WHERE orgId = ? AND expGroup = ? AND condition = ?",
                                          {}, $orgId, $exp->{expGroup}, $exp->{condition_1});
    foreach my $row (@$specOG) {
        my ($locusId,$nInOG) = @$row;
        $cons{$locusId} = $nInOG;
    }
}

print $cgi->h3($header) if defined $header;

if (@fit > 0) { # show the table
  my @trows = ();
  foreach my $row (@fit) {
    my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$row->{locusId}&expName=$expName";
    $strainUrl .= "&help=1" if $help;
    my $orthUrl = "orthFit.cgi?orgId=$orgId&locusId=$row->{locusId}"
      . "&expGroup=" . uri_escape($exp->{expGroup})
        . "&condition1=" . uri_escape($exp->{condition_1});
    $orthUrl .= "&help=1" if $help;
    push @trows,
      $cgi->Tr( {align => 'left', valign => 'top'},
                td(checkbox('locusId',0,$row->{locusId},'')),
                td(Utils::gene_link($dbh, $row, "name", "myFitShow.cgi")),
                td($row->{gene}),
                td({ -bgcolor => Utils::fitcolor($row->{fit}) },
                   a({ -href => $strainUrl,
                       -title => "per-strain data",
                       -style => "color:rgb(0,0,0)" },
                     sprintf("%.1f", $row->{fit})) ),
                td( sprintf("%.1f", $row->{t}) ),
                td(Utils::gene_link($dbh, $row, "desc", "domains.cgi")),
                td(a({ title => "Compare to data from similar experiments or orthologs",
                       href => $orthUrl},
                     exists $cons{$row->{locusId}} && $cons{$row->{locusId}} > 1?
                     "<i>conserved</i>" : "compare")) );
  }
  print
    start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
      hidden('orgId', $orgId),
	$cgi->table( { cellspacing => 0, cellpadding => 3},
                     $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                              $cgi->th(['&nbsp;', 'gene','name','fitness','t score','description', '&nbsp;']),
                              @trows) ),
                                "<BR><BR>",
                                  submit(-name => 'Heatmap of selected genes'),
                                    end_form;
} elsif ($show eq "quality") {
    my @trows = ();
    push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			  $cgi->td([ "Time0", $exp->{timeZeroSet}, "which Time0s the sample was compared to"]));
    my %expl = ("cor12" => "rank correlation(fit1, fit2), where fit1 is fitness for the first half (10-50%) and fit2 is fitness for the second half (50-90%) of each gene",
		"maxFit" => "The maximum fitness value",
		"opcor" => "rank correlation(upstream gene, downstream gene) over pairs that are adjacent and likely to be in the same operon",
		"adjcor" => "like opcor but for adjacent genes that are not on the same strand",
		"gccor" => "linear correlation of gene fitness and gene GC content",
		"mad12" => "median absolute difference of fit1, fit2",
		"mad12c" => "median absolute difference of log count for 1st and 2nd half of genes in this sample",
		"mad12c_t0" => "like mad12c but for the Time0s",
		"gMed" => "median reads per gene in this sample",
		"gMedt0" => "median reads per gene in the Time0 sample",
		"gMean" => "mean reads per gene in this sample",
		"nMapped" => "#reads for this sample that corresponded to a known strain (in millions)",
		"nPastEnd" => "#reads that corresponded to a strain that has an insertion within the suicide vector instead of within the genome.",
		"nGenic" => "#reads that lie within central 10-90% of a gene",
		"nUsed" => "#reads used to estimate gene fitness (genic and enough coverage for strain and for gene)" );
    foreach my $key (qw{cor12 maxFit opcor adjcor gccor mad12 mad12c mad12c_t0}) {
	my $val = $exp->{$key};
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.2f", $exp->{$key}), $expl{$key}]));
    }
    foreach my $key (qw{gMed gMedt0 gMean}) {
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.0f", $exp->{$key}), $expl{$key}]));
    }
    foreach my $key (qw{nMapped nPastEnd nGenic nUsed}) {
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.3f M", $exp->{$key}/1e6), $expl{$key}]));
    }
    print $cgi->table({cellpadding => 3, cellspacing => 0}, @trows);
}
print "\n";

if ($show ne "specific") {

    print h3("Specific Phenotypes");
    if (@$spec > 0) {
	print p(a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=specific"},
                  "For " . scalar(@$spec). " genes in this experiment"));
        print p(a({href => "spec.cgi?orgId=$orgId&expGroup=".uri_escape($exp->{expGroup})."#".$exp->{condition_1} },
                  "For $exp->{expGroup} $exp->{condition_1} in $orginfo->{$orgId}{genome}"))
            if $exp->{expGroup} && $exp->{condition_1};
    }  else {
	print $cgi->p("None in this experiment");
        if ($exp->{expGroup}) {
            print
                p(a({href => "spec.cgi?orgId=$orgId&expGroup=" . uri_escape($exp->{expGroup})
                         . ($exp->{condition_1} eq "" ? "" : "#" . uri_escape($exp->{condition_1}))},
                    "For $orginfo->{$orgId}{genome} in $exp->{expGroup} experiments"));
        }
    }
}
print "\n";

if ($exp->{condition_1} ne "") {
    print
        p(a({href => "orthCond.cgi?expGroup=" . uri_escape($exp->{expGroup})
                 . "&condition1=" . uri_escape($exp->{condition_1})},
            ($show eq "specific" ? "Specific phenotypes for" : "For")
            . " $exp->{expGroup} $exp->{condition_1} across organisms"));
}
print "\n";

if ($show eq "") {

  if (@$spec > 0) {
    # Show the relevant SEED subsystems
    # Need to use SEEDAnnotation => SEEDAnnotationToRoles => SEEDRoles.subsystem
    my %subsys = (); # seed subsystem => # of genes
    foreach my $gene (@$spec) {
      my $subsys = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT subsystem FROM SeedAnnotation
                                                JOIN SeedAnnotationToRoles USING (seed_desc)
                                                JOIN SeedRoles USING (seedrole)
                                                WHERE orgId = ? AND locusId = ? },
                                            {}, $orgId, $gene->{locusId});
      foreach my $subsys (@$subsys) {
        $subsys{$subsys}++;
      }
    }
    my @subsys = sort { $subsys{$b} <=> $subsys{$a} || $a cmp $b } keys %subsys;
    if (@subsys > 0) {
      print h3("SEED Subsystems");
      my @th = (th("Subsystem"), th("#Specific"));
      my @trows = ( Tr(@th) );
      foreach my $subsys (@subsys) {
        my $subsysShow = $subsys; $subsysShow =~ s/_/ /g;
        my @trow = (td(a({ -href => "seedsubsystem.cgi?orgId=$orgId&subsystem=$subsys&expName=$expName" }, $subsysShow)),
                    td($subsys{$subsys}));
        push @trows, Tr(@trow);
      }
      print table( { -cellpadding => 3, -cellspacing => 0}, @trows);
    }

    # try to highlight useful maps that contain genes with specific phenotypes
    my %ec = ();
    my %rxn = ();
    foreach my $gene (@$spec) {
      my $locusId = $gene->{locusId};
      $ec{$locusId} = Utils::GeneToEc($dbh, $orgId, $locusId);
      $rxn{$locusId} = Utils::GeneToRxn($dbh, $orgId, $locusId, $ec{$locusId});
    }
    my %specEc = (); # ec => 1 for those locusIds
    foreach my $list (values %ec) {
      foreach my $ec (@$list) {
        $specEc{$ec} = 1;
      }
    }
    my %specRxn = ();
    foreach my $list (values %rxn) {
      foreach my $rxn (@$list) {
        $specRxn{$rxn} = 1;
      }
    }

    print h3("Metabolic Maps"),
      p("Color code by fitness: see",
        a({href => "keggmap.cgi?mapId=01100&orgId=$orgId&expName=$expName"}, "overview map"),
        "or",
        a({href => "keggmaplist.cgi?orgId=$orgId&expName=$expName"}, "list of maps")."." );
    if (keys(%specEc) > 0) {
      my $ecIn = join(",", map "'".$_."'", keys %specEc);
      my $maps = $dbh->selectall_arrayref(qq{ SELECT mapId, title, COUNT(DISTINCT objectId) nEc
                                                FROM KEGGConf JOIN KEGGMap USING (mapId)
                                                WHERE objectId IN ( $ecIn ) AND type=1
                                                GROUP BY mapId
                                                ORDER BY nEc DESC });
      if (scalar(@$maps) > 0) {
        my @mapShow = ();
        foreach my $map (@$maps) {
          my ($mapId,$title,$nEc) = @$map;
          push @mapShow, li(a({href => "keggmap.cgi?orgId=$orgId&expName=$expName&mapId=$mapId"},
                              $title));
        }
        print p("Maps containing gene(s) with specific phenotypes:",
                  ul(@mapShow));
      }
    }
    print "\n";

    # Map reactions to pathways
    if (keys %specRxn > 0) {
      my $rxnSpec = join(",", map { "'" . $_ . "'" } keys %specRxn);
      my $rxnPath = $dbh->selectall_arrayref(qq{ SELECT * FROM MetaCycPathwayReaction
                                                 JOIN MetaCycPathway USING (pathwayId)
                                                 WHERE rxnId IN ($rxnSpec) },
                                             { Slice => {} });
      my %path = (); # pathway to hash of pathwayId, pathwayName, nSpecific, etc.
      foreach my $row (@$rxnPath) {
        my $pathId = $row->{pathwayId};
        unless (exists $path{$pathId}) {
          $path{$pathId} = { "pathwayId" => $pathId, "pathwayName" => $row->{pathwayName}, "nSpecific" => 0 };
        }
        $path{$pathId}{nSpecific}++;
      }
      my @path = values %path;
      if (@path > 0) {
        foreach my $row (@path) {
          my $cov = $dbh->selectrow_hashref("SELECT * FROM MetacycPathwayCoverage WHERE orgId = ? AND pathwayId = ?",
                                             { Slice => {} }, $orgId, $row->{pathwayId});
          foreach my $field (qw{nSteps nFound}) {
            $row->{$field} = $cov->{$field};
          }
        }

        my @trows = ();
        push @trows, $cgi->Tr(th(["Pathway", "#Steps", "#Present", "#Specific"]));
        @path = sort { $b->{nSpecific} / $b->{nSteps} <=> $a->{nSpecific} / $a->{nSteps}
                         || $b->{nFound} / $b->{nSteps} <=> $a->{nFound} / $a->{nSteps}
                           || $b->{nSteps} <=> $a->{nSteps} } @path;
        foreach my $row (@path) {
          push @trows,
            $cgi->Tr(td([ a({ -href => "pathway.cgi?orgId=$orgId&expName=$expName&pathwayId=$row->{pathwayId}" },
                            $row->{pathwayName}),
                          $row->{nSteps},
                          $row->{nFound},
                          $row->{nSpecific} ]));
        }
        print
          h3("MetaCyc Pathways"),
            p(a({ -href => "pathwaysOrg.cgi?orgId=$orgId" }, "Pathways"),
              "that contain genes with specific phenotypes:"),
              $cgi->table({cellpadding => 3, cellspacing => 0}, @trows),
                "\n";
      }
    }
  }
}
$dbh->disconnect();
Utils::endHtml($cgi);

sub CompoundToHTML($$$$) {
    my ($compound, $concentration, $units, $mediaOrMix) = @_;
    return "" if $compound eq "";
    return "$concentration$units $compound" if exists $mediaOrMix->{$compound};
    my ($cas) = $dbh->selectrow_array("SELECT CAS FROM Compounds WHERE compound=?", {}, $compound);
    $compound =~ s/^supernat[ea]nt; /supernatant from /i;
    my $html = $cas ? a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$cas"}, $compound)
      : $compound;
    return $html if $concentration eq "" && $units eq "";
    return "$html ($concentration $units)";
}

sub MediaComponentHTML($$) {
  my ($components, $strength) = @_;
  # mix => list of components; note mix = "" for most
  my %compStrings = ("" => []);
  foreach my $row (@$components) {
    my $compString = $row->{CAS} ?
      a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$row->{CAS}"}, $row->{compound})
        :  $row->{compound};
    if ($row->{concentration} && $row->{units}) {
      my $conc = $row->{concentration} * $strength;
      $compString = join(" ", $conc, $row->{units}, $compString);
      $row->{mix} = "" if !defined $row->{mix};
      push @{ $compStrings{$row->{mix}} }, $compString;
    }
  }
  # The non-mix components
  my $html = join(", ", @{ $compStrings{""} });
  # The mix components get a label
  foreach my $mix (sort keys %compStrings) {
    next if $mix eq "";
    $html .= ", $mix " . small("(" . join(", ", @{ $compStrings{$mix} }) . ")");
  }
  return $html;
}
