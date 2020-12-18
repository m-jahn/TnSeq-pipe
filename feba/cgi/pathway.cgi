#!/usr/bin/perl -w
#######################################################
## pathway.cgi -- view MetaCyc pathway with candidate
##	genes in that organism and fitness data
##
## Copyright (c) 2018 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# pathwayId -- which MetaCyc pathway
#
# Optional CGI parameters:
# expName -- which experiment(s) to show fitness data from
# addexp -- add experiment(s) which match this query

use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $pathId = $cgi->param('pathwayId') || die "No pathwayId parameter";
my @expNames = $cgi->param('expName');
my $addexp = $cgi->param('addexp');
my $debug = $cgi->param('debug') || 0;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

my ($pathwayName) = $dbh->selectrow_array("SELECT pathwayName FROM MetacycPathway WHERE pathwayId = ?",
                                          {}, $pathId);
die "Unknown pathway: $pathId" unless defined $pathwayName;

# Fetch metadata on experiments
my %exps = map { $_->{expName} => $_ } @{ $dbh->selectall_arrayref("SELECT * FROM Experiment WHERE orgId = ?",
                                                                   { Slice => {} }, $orgId) };

my @warnings = ();
if ($addexp) {
  my %expSeen = map { $_ => 1 } @expNames;
  my $addexps = Utils::matching_exps($dbh, $orgId, $addexp);
  foreach my $exp (@$addexps) {
    my $expName = $exp->{expName};
    next if exists $expSeen{$expName};
    $expSeen{$expName} = 1;
    push @expNames, $expName;
  }
  push @warnings, span({ -style => "color: red;" },
                       qq{No experiments matching "$addexp"})
    if @$addexps == 0;
}

my @exps = ();
foreach my $expName (@expNames) {
  die "Invalid expName $expName" unless $expName =~ m/^[a-zA-Z][a-zA-Z0-9_-]+$/;
  die "Invalid expName $expName for org $orgId" unless exists $exps{$expName};
  push @exps, $exps{$expName};
}

print
  header,
  Utils::start_page("MetaCyc Pathway: $pathwayName in $orginfo->{$orgId}{genome}"),
  h2("MetaCyc Pathway: ",
     a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=$pathId" }, $pathwayName),
     "in",
     a({ -href => "org.cgi?orgId=$orgId" }, $orginfo->{$orgId}{genome})),
  p(join(br(), @warnings)),
  # the experiment selector
  start_form(-name => 'input', -method => 'GET', -action => 'pathway.cgi'),
  hidden( -name => 'orgId', -value => $orgId, -override => 1),
  hidden( -name => 'pathwayId', -value => $pathId, -override => 1),
  join("\n", map { hidden( -name => 'expName', -value => $_, -override => 1) } @expNames),
  p({-class => "buttons", style=>"align:left; white-space:nowrap; line-height:40px;"}, "Add experiment(s): ",
    textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
    submit('Add','Add') ),
  end_form,
  br();

my $rxns = $dbh->selectall_arrayref(qq{ SELECT * FROM MetacycPathwayReaction
                                        JOIN MetacycReaction USING (rxnId)
                                        WHERE pathwayId = ? },
                                    { Slice => {} }, $pathId);
die "Empty pathway: $pathId\n" unless @$rxns > 0;
my %rxns = map { $_->{rxnId} => $_ } @$rxns;

my $preds = $dbh->selectall_arrayref("SELECT * FROM MetacycPathwayReactionPredecessor WHERE pathwayId = ?",
                                    { Slice => {} }, $pathId);
# %pred has rxnId => list of predecessor rxnIds, only including those in this pathway
my %pred = map { $_ => [] } keys %rxns;
foreach my $row (@$preds) {
  push @{ $pred{$row->{rxnId}} }, $row->{predecessorId}
    if exists $rxns{ $row->{predecessorId} } && exists $rxns{ $row->{rxnId} };
}

# Order the reactions by precedence. (hack-ish)

# First, every node with no precedors gets 0
my %score = ();
while (my ($rxnId, $predlist) = each %pred) {
  if (scalar(@$predlist) == 0) {
    $score{$rxnId} = 0;
  }
}

# For strictly circular pathways such as GLYOXYLATE-BYPASS, this will fail to choose a starting point, so, give some node a score
if (scalar(keys %score) == 0) {
  my $first = (sort keys %pred)[0];
  $score{$first} = 0;
}

# Then, for every node with all precedors having a score, set to max(preceding score) + 1
# If this fails to make any assignments, and not all scores are set, then relax this to allow just 1 score to be set
#	(This is a way to get around cycles)
my $nRound = 0;
my $maxRounds = 100;
my $strict = 1;
for($nRound = 0; $nRound < $maxRounds; $nRound++) {
  my $nSet = 0;
  while (my ($succId, $predlist) = each %pred) {
    next if exists $score{$succId}; # never reassign a score
    my $nPredSet = 0;
    my $maxPred = 0;
    foreach my $rxnId (@$predlist) {
      if (exists $score{$rxnId}) {
        $nPredSet++;
        $maxPred = $score{$rxnId} if $maxPred < $score{$rxnId};
      }
    }
    if ($nPredSet == scalar(@$predlist) || ($nPredSet > 0 && ! $strict)) {
      $score{$succId} = $maxPred + 1;
      $nSet++;
    }
  }
  if ($nSet == 0) {
    last unless $strict;
    $strict = 0;
    print p("Strict off") if $debug;
  } else {
    $strict = 1;
    print p("Strict on after setting $nSet") if $debug;
  }
}
print p("Sorting the pathway failed to terminate") if $nRound == $maxRounds;

# some reactions may still not have scores if there are subgraph cycles; those
# will be placed arbitrarily
my @rxnIds = sort { ($score{$a} || 0) <=> ($score{$b} || 0) || $a cmp $b } keys %rxns;

if ($debug) {
  print p("Reaction ordering and scores after $nRound rounds:");
  foreach my $rxnId (@rxnIds) {
    print p($rxnId, "score", $score{$rxnId}, "predecessors", @{ $pred{$rxnId} });
  }
}

my $primary = $dbh->selectall_arrayref("SELECT * FROM MetacycPathwayPrimaryCompound WHERE pathwayId = ?",
                                       { Slice => {} }, $pathId);
my %primary = (); # rxnId => compoundId => side => 1
foreach my $row (@$primary) {
  $primary{ $row->{rxnId} }{ $row->{compoundId} }{ $row->{side} } = 1;
}

my @trows = ();
if (@exps > 0) {
  # make a header row
  my @hrow = ( th("Reactions and Genes") );
  foreach my $exp (@exps) {
    push @hrow, th(a( { -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}",
                       -title => $exp->{expName} },
                     $exp->{expDesc} ));
  }
  push @trows, Tr(@hrow)."\n";
  # make a removal row
  my @rmrow = ( td("") );
  my $base = "pathway.cgi?orgId=$orgId&pathwayId=$pathId";
  foreach my $expId (@expNames) {
    my @expNames2 = grep { $_ ne $expId } @expNames;
    my $URL = join("&", $base, map { "expName=$_" } @expNames2);
    push @rmrow, td(a( { -href => $URL}, "remove"));
  }
  push @trows, Tr(@rmrow)."\n";
}

my $genestyle = "padding-left: 3em;"; # inset gene descriptions

my %primaryShown = (); # compoundId => 1 if shown as primary already
my $ncol = 1 + scalar(@expNames);
# rxnId => {isSpontaneous => 0 or 1, loci => list of locusIds }
my $cands = Utils::MetacycPathwayCandidates($dbh, $orgId, $pathId);
my @loci = (); # list of loci shown
foreach my $rxnId (@rxnIds) {
  my $rxn = $rxns{$rxnId};
  my $loci = $cands->{$rxnId}{loci};
  my $ecs = $dbh->selectcol_arrayref("SELECT ecnum FROM MetacycReactionEC WHERE rxnId = ?",
                                    {}, $rxnId);
  my $ecdesc = "";
  ($ecdesc) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ? LIMIT 1",
                                     {}, $ecs->[0])
    if @$ecs > 0;
  my @eclinks = map a({-href => "myFitShow.cgi?orgId=$orgId&gene=ec:$_" }, $_), @$ecs;

  # Ignore reaction names that are just EC numbers
  my $rxnName = $rxn->{rxnName};
  $rxnName = "" if $rxnName =~ m/^[0-9][.][0-9]+[.]/;
  $rxnName = $ecdesc if $rxnName eq "" && defined $ecdesc;
  $rxnName =~ s/[.]$//;
  $rxnName .= " (in reverse)" if $rxnName && $rxn->{direction} == -1;

  my @left = ();
  my @right = ();
  my $cmps = $dbh->selectall_arrayref(qq{ SELECT * FROM MetacycReaction
                                         JOIN MetacycReactionCompound USING (rxnId)
                                         LEFT JOIN MetacycCompound USING (compoundId)
                                         WHERE rxnId = ? },
                                     { Slice => {} }, $rxnId);
  foreach my $i (0..(scalar(@$cmps)-1)) {
    $cmps->[$i]{row} = $i;
  }
  my @sortedcmps = sort { (exists $primary{$rxnId}{$b->{compoundId}}) - (exists $primary{$rxnId}{$a->{compoundId}})
                            || (exists $primaryShown{$b->{compoundId}}) - (exists $primaryShown{$a->{compoundId}})
                              || $a->{row} - $b->{row} } @$cmps;
  foreach my $row (@sortedcmps) {
    my $compoundId = $row->{"compoundId"};
    # Some reactants are actually compound classes -- for those, just use the name.
    my $name = $row->{compoundName} || $compoundId;
    my $color = undef;
    if (exists $primary{$rxnId}{$compoundId}) {
      $primaryShown{$compoundId} = 1;
    } else {
      $color = "DimGrey";
    }
    $name = span({ style => "color: $color" }, $name) if defined $color;
    my $side = $rxn->{direction} * $row->{"side"};
    my $list = $side == -1 ? \@left : \@right;
    $name = $row->{coefficient} . " " . $name unless $row->{coefficient} eq "" || $row->{coefficient} eq "1";
    $name = $name . Sub("[" . lc($row->{compartment}) . "]") if $row->{compartment} ne "";
    push @$list, $name;
  }
  my $spontaneous = $rxn->{isSpontaneous} ? "(spontaneous)" : "";

  my @locirows = ();
  # show list of genes
  foreach my $locusId (@$loci) {
    my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
                                       {}, $orgId, $locusId);
    my @generow = ();
    push @generow, td({ -style => $genestyle },
                      Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi")
                     . ($gene->{gene} ? " (" . $gene->{gene} . ")" : ""));
    push @loci, $gene->{locusId};
    my $fitrows = $dbh->selectall_arrayref("SELECT expName,fit,t FROM GeneFitness WHERE orgId=? AND locusId=?",
                                       {}, $orgId, $locusId);
    my %fit = map { $_->[0] => $_ } @$fitrows;
    foreach my $expName (@expNames) {
      my $fitrow = $fit{$expName};
      my ($fit,$t);
      (undef,$fit,$t) = @$fitrow if defined $fitrow;
      push @generow, td({ -bgcolor => Utils::fitcolor($fit) },
                        a({ -href => "strainTable.cgi?orgId=$orgId&locusId=$locusId&expName=$expName",
                            -title => defined $fitrow ? sprintf("t = %.1f", $t) : "No data",
                            -style => "color:rgb(0,0,0)" },
                          defined $fitrow ? sprintf("%.1f", $fit) : ""));
    }
    push @locirows, Tr(@generow);
  }
  push @trows, Tr(td({ -colspan => $ncol},
                     a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=REACTION&object=$rxnId",
                         -title => "MetaCyc Reaction" },
                       ($rxnName ? $rxnName . ":" . br() : "")
                       . join(" + ", @left) . "&rarr;" . join(" + ", @right)),
                     $spontaneous,
                     @$ecs > 0 ? "(EC " . join("; ", @eclinks) . ")" : ""));
  if (@locirows > 0) {
    push @trows, @locirows;
  } else {
    push @trows, Tr(td({ -style => $genestyle, -colspan => $ncol }, "No genes"));
  }
}
print "\n", table( { cellspacing => 0, cellpadding => 3 }, @trows);

print p("Links:"),
  start_ul(),
  li( a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=$pathId&detail-level=3" }, "Pathway details at MetaCyc"), );

my @lociUniq = ();
my %lociSeen = ();
foreach my $locusId (@loci) {
  push @lociUniq, $locusId unless exists $lociSeen{$locusId};
  $lociSeen{$locusId} = 1;
}
# When written, this should always be the case because pathways with no candidate
# genes should not be linked to, but just in case
if (@lociUniq > 0) {
  my $URL = "genesFit.cgi?orgId=$orgId&showAll=0&" . join("&", map { "locusId=$_" } @lociUniq);
  my $nLoci = scalar(@lociUniq);
  print li( a({ -href => $URL }, "Fitness data for $nLoci genes") );
}

print
  li( start_form(-name => 'orgselect', -method => 'GET', -action => 'pathway.cgi'),
     hidden( -name => 'pathwayId', -value => $pathId, -override => 1),
     "Select organism:",
     Utils::OrgSelector($orgId, $orginfo),
     "<button type='submit'>Go</button>",
     end_form ),
  li( a({ -href => "pathwaysOrg.cgi?orgId=$orgId" }, "All MetaCyc pathways for $orginfo->{$orgId}{genome}") ),
  end_ul(),
  p(small(qq{Candidate genes for each reaction are identified from best hits to MetaCyc; by matching EC numbers (which are assigned by TIGRFam, SEED, or best hits to KEGG); or by matching SEED roles to KEGG reactions to MetaCyc reactions. See the "Protein" tab of each gene for more information}));


$dbh->disconnect();
Utils::endHtml($cgi);

