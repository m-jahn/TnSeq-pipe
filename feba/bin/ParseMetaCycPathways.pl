#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use ParsePTools;

# Parse pathways, reactions, and compounds from MetaCyc's data directory

my $out = ".";

my $usage = <<END
Usage: ParseMetaCycPathways.pl [ -out $out ] -in metacyc_data_directory

Parses flat files in the MetaCyc data directory to produce
tab-delimited files in the output directory:
MetacycPathway.tab
MetacycPathwayReaction.tab
MetacycPathwayReactionPredecessor.tab
MetacycPathwayPrimaryCompound.tab
MetacycReaction.tab
MetacycReactionCompound.tab
MetacycCompound.tab

These have no header line and should match the schema, see
feba/lib/db_setup_tables.sql
END
;

my $indir;
die $usage unless GetOptions('out=s' => \$out,
                             'in=s' => \$indir)
  && @ARGV == 0;
die "Must specify -in:\n$usage" unless defined $indir;
foreach my $dir ($indir, $out) {
  die "Not a directory: $dir\n" unless -d  $dir;
}

foreach my $stub (qw/pathways reactions compounds/) {
  die "No such file: $indir/$stub.dat\n"
    unless -e "$indir/$stub.dat";
}

my %compounds = (); # compoundId => hash of compoundName, keggLigand, formula (as string)
open (my $fhc, "<", "$indir/compounds.dat")
  || die "Cannot open $indir/compounds.dat";
while (my $cmp = ParsePTools($fhc)) {
  my $cmpId = $cmp->{"UNIQUE-ID"}[0]{"value"};
  die "Invalid compound in $indir/compounds.dat" unless $cmpId;
  my $name = $cmp->{"COMMON-NAME"}[0]{"value"} || "";
  my @formula = ();
  foreach my $l (@{ $cmp->{"CHEMICAL-FORMULA"} }) {
    my $string = $l->{"value"};
    die "Invalid component of chemical formula $string" unless $string =~ m/^[(]([A-Za-z-]+) (\d+)[)]$/;
    my ($element, $cnt) = ($1,$2);
    push @formula, ucfirst(lc($element)) . $cnt;
  }
  my $keggLigand = "";
  foreach my $l (@{ $cmp->{"DBLINKS"} }) {
    my $string = $l->{"value"};
    $keggLigand = $1 if $string =~ m/^[(]LIGAND-CPD "([A-Z0-9]+)"/;
  }
  die "Duplicate compound $cmpId" if exists $compounds{$cmpId};
  $compounds{$cmpId} = { "compoundName" => $name,
                         "keggLigand" => $keggLigand,
                         "formula" => join("", @formula) };
}
close($fhc) || die "Error reading $indir/compounds.dat";
print STDERR "Read " . scalar(keys %compounds) . " compounds\n";

# Classes are hierarchical, so a compound class is not directly labelled as such --
# one needs to infer that it is a compound by going up the hierarchy
# (i.e. All-Carbohydrates has TYPES = Compounds and Glycans has TYPES = Carbohydrates.)
# Instead, will just assume that it is a compound if a pathway or reaction refers to it
my %classes = (); # id => hash that includes id, name and types (also a hash)
open(my $fhcc, "<", "$indir/classes.dat")
  || die "Cannot open $indir/classes.dat";
while(my $cc = ParsePTools($fhcc)) {
  my $id = $cc->{"UNIQUE-ID"}[0]{"value"};
  die unless $id;
  die "Duplicate class $id" if exists $classes{$id};
  my %obj = ( "id" => $id, "types" => {}, "name" => "" );
  foreach my $l (@{ $cc->{"TYPES"} }) {
    $obj{"types"}{ $l->{"value"} } = 1;
  }
  $obj{"name"} = $cc->{"COMMON-NAME"}[0]{"value"}
    if exists $cc->{"COMMON-NAME"};
  $classes{$id} = \%obj;
}
close($fhcc) || die "Error reading $indir/classes.dat";

my %path = (); # pathwayId => hash of pathwayName, subpathways, reactions, hypo_reactions, and predecessor
# subpathways is a list of pathway ids
# reactions is a hash of rxnId => hash with fields rxnId, direction, primary,
#	where primary is a list of [compoundId,side]
#	subpathways may also end up in the reaction list but these will be ignored later on.
# hypo_reactions maps rxnId => 1
# predecessor is a hash of rxnId to list of predecessor rxnIds
my %dirToVal = ("NIL" => "", ":L2R" => 1, ":R2L" => -1);
my %otherCompounds = ();
open (my $fhpath, "<", "$indir/pathways.dat")
  || die "Cannot open $indir/pathways.dat";
while (my $path = ParsePTools($fhpath)) {
  my $pathwayId = $path->{"UNIQUE-ID"}[0]{"value"};
  die unless $pathwayId;
  die "Duplicate pathway $pathwayId" if exists $path{$pathwayId};
  my $pathwayName = $path->{"COMMON-NAME"}[0]{"value"} || "";
  #my %types = map { $_->{"value"} => 1 } @{ $path->{"TYPES"} };
  my %pred = ();

  # Parse REACTION-LAYOUT to get the reactions, then set up isHypothetical and predecessor(s)
  my %reactions = ();
  foreach my $l (@{ $path->{"REACTION-LAYOUT"} }) {
    my $string = $l->{"value"};
    $string =~ m/^[(]([^() ]+) +[(]:LEFT-PRIMARIES *([^()]*)[)] +[(]:DIRECTION ([^ ()]+)[)] [(]:RIGHT-PRIMARIES *([^()]*)[)][)]$/
      || die "Cannot parse REACTION-LAYOUT $string";
    my ($rxnId, $lefts, $dirstring, $rights) = ($1,$2,$3,$4);
    my $dir = $dirToVal{$dirstring};
    die "Unknown direction $dir for pathway $pathwayId REACTION-LAYOUT $string"
      unless defined $dir;
    my @lefts = split / +/, $lefts;
    my @rights = split / +/, $rights;
    # These should be present in the compounds table or the classes table
    my @both = @lefts;
    push @both, @rights;
    foreach my $id (@both) {
      $otherCompounds{$id} = 1 unless exists $compounds{$id};
    }
    my @primary = map [$_, -1], @lefts;
    push @primary, map [$_, +1], @rights;
    $reactions{$rxnId} = { "rxnId" => $rxnId,
                           "direction" => $dir,
                           "primary" => \@primary };
  }

  foreach my $l (@{ $path->{"PREDECESSORS"} }) {
    my $string = $l->{"value"};
    # a pathway can be listed as a predecessor; this means incorporate all predecessors from that pathway,
    # but will rely on SUB-PATHWAYS to explicitly incorporate subpathways
    if ($string =~ m/^[A-Z0-9+-]+$/) {
      # ignore;
    } else {
      $string =~ m/^[(]([^()]+)[)]$/ || die "Cannot parse PREDECESSORS for pathway $pathwayId: $string";
      my $pre = $1;
      # both reaction and predecessor ids are often quoted, not sure why
      $pre =~ s/"//g;
      my @pre = split / +/, $pre;
      my $rxnId = shift @pre;
      push @{ $pred{$rxnId} }, @pre;
    }
  }

  my %hypo_reactions = ();
  foreach my $l (@{ $path->{"HYPOTHETICAL-REACTIONS"} }) {
    my $rxnId = $l->{"value"};
    $hypo_reactions{$rxnId} = 1;
  }

  my %subpathways = ();
  foreach my $l (@{ $path->{"SUB-PATHWAYS"} }) {
    $subpathways{ $l->{value} } = 1;
  }
  my @subpathways = sort keys %subpathways;

  my @reactions = sort { $a->{"rxnId"} cmp $b->{"rxnId"} } values(%reactions);
  $path{$pathwayId} = { "pathwayName" => $pathwayName,
                        "reactions" => \@reactions,
                        "subpathways" => \@subpathways,
                        "hypo_reactions" => \%hypo_reactions,
                        "predecessor" => \%pred };
}
close($fhpath) || die "Error reading $indir/pathways.dat";
print STDERR "Read " . scalar(keys %path) . " base-level pathways with "
  . scalar(keys %otherCompounds) . " unknown compounds (classes?)\n";
foreach my $id (keys %otherCompounds) { # this is rare
  print STDERR "Warning: unknown compound $id is not a class\n"
    unless exists $classes{$id};
}

# Recursively expand subpathways
my %path_expanded = (); # pathway => subpathway => 1 if already incorporated
for(;;) {
  my $nChange = 0;
  while (my ($pathId, $path) = each %path) {
    foreach my $subId (@{ $path->{subpathways} }) {
      next if exists $path_expanded{$pathId}{$subId};
      # In principle, only actual pathways should end up in this list,
      # but this is not always followed;
      next unless exists $path{$subId};
      my $sub = $path{$subId};
      my $grandIds = $sub->{subpathways};
      # only expand sub into path if sub is fully expanded already
      my $fully = 1;
      foreach my $grandId (@$grandIds) {
        $fully = 0 if exists $path{$grandId} && !exists $path_expanded{$subId}{$grandId};
      }
      if ($fully) {
        # expand sub into path:
        $path_expanded{$pathId}{$subId} = 1;
        $nChange++;
        # Add all the reactions
        push @{ $path->{reactions} }, @{ $sub->{reactions} };
        # Add all the predecessor relationships
        while (my ($rxnId, $pre) = each %{ $sub->{predecessor} }) {
          push @{ $path->{predecessor}{$rxnId} }, @$pre;
        }
        # Add isHypothetical (not sure if this is necessary, it may be represented in the parent already)
        foreach my $rxnId (keys %{ $sub->{hypo_reactions} }) {
          $path->{hypo_reactions}{$rxnId} = 1;
        }
      }
    }
  }
  print STDERR "Expanded $nChange subpathways\n" if $nChange > 0;
  last if $nChange == 0;
}

# rxnId => hash of rxnId, rxnName, ecnum (a list), isSpontaneous, keggrxnId, compounds, subreactions (a list of ids)
#     compounds is a list of hashes of compoundId, side, coefficient, compartment
# Note that a few reactions are compound reactions, as indicated using REACTION-LIST and stored in subreactions
my %rxn = ();

open (my $fhr
, "<", "$indir/reactions.dat")
  || die "Cannot open $indir/reactions.dat";
while (my $rxn = ParsePTools($fhr)) {
  my $rxnId = $rxn->{"UNIQUE-ID"}[0]{"value"};
  die "Missing UNIQUE-ID" unless $rxnId;
  die "Duplicate reaction id $rxnId" if exists $rxn{$rxnId};

  my $obj = { "rxnId" => $rxnId, "ecnum" => [], "isSpontaneous" => 0, "compounds" => [],
              "keggrxnId" => "", "subreactions" => [] };

  foreach my $l (@{ $rxn->{"EC-NUMBER"} }) {
    my $ec = $l->{"value"};
    # sometimes EC number is like |EC-1.8.3.b| instead of EC-1.8.3.b
    $ec =~ s/^[|]//;
    $ec =~ s/[|]$//;
    $ec =~ s/^EC-//;
    my @ec = split /[.]/, $ec;
    while (@ec < 4) {
      push @ec, '-';
    }
    push @{ $obj->{"ecnum"} }, join(".", @ec);
  }

  if (exists $rxn->{"SPONTANEOUS?"} && $rxn->{"SPONTANEOUS?"}[0]{"value"} eq "T") {
    $obj->{"isSpontaneous"} = 1;
  }

  foreach my $l (@{ $rxn->{"REACTION-LIST"} }) {
    push @{ $obj->{subreactions} }, $l->{value};
  }

  my @cmp = ();
  foreach my $l (@{ $rxn->{"LEFT"} }) {
    $l->{"side"} = -1;
    push @cmp, $l;
  }
  foreach my $l (@{ $rxn->{"RIGHT"} }) {
    $l->{"side"} = +1;
    push @cmp, $l;
  }

  # ensure that a compound only shows up once on each side/compartment
  # first transform to the list into the final object form
  my @cmp2 = ();
  foreach my $l (@cmp) {
    my $cmpId = $l->{"value"};
    my $coefficient = $l->{"COEFFICIENT"} || 1;
    my $compartment = $l->{"COMPARTMENT"} || "";
    $compartment =~ s/^CCO-//;
    push @cmp2, { "compoundId" => $cmpId,
                  "side" => $l->{"side"},
                  "coefficient" => $coefficient,
                  "compartment" => $compartment };
    $otherCompounds{$cmpId} = 1 unless exists $compounds{$cmpId};
  }

  # then identify duplicate compounds
  my %rcmp = ();
  foreach my $c (@cmp2) {
    my $key = join(":::", $c->{"side"}, $c->{"compoundId"}, $c->{"compartment"});
    push @{ $rcmp{$key} }, $c;
  }
  # and sum coefficient (if possible)
  foreach my $key (sort keys %rcmp) {
    my @k = @{ $rcmp{$key} };
    my $l = $k[0];
    if (@k > 1) {
      # try to combine the coefficients (all other attributes match)
      my @coeff = map $_->{"coefficient"}, @k;
      # are they all numbers?
      my @numcoeff = grep m/^-?\d+[.]?\d*$/, @coeff;
      if (scalar(@coeff) == scalar(@numcoeff)) {
        my $tot = 0;
        foreach my $coeff (@coeff) {
          $tot += $coeff;
        }
        $l->{"coefficient"} = $tot;
      } else {
        # just use the first value
        print STDERR "Warning: duplicate compound " . $key->{"compoundId"} . " without numeric coefficients in reaction $rxnId\n";
      }
    }
    push @{ $obj->{"compounds"} }, $l;
  }

  # and identify a link to a KEGG reaction, if any
  foreach my $l (@{ $rxn->{"DBLINKS"} }) {
    my $link = $l->{"value"};
    if ($link =~ m/^[(]LIGAND-RXN "([A-Z0-9]+)"/) {
      $obj->{keggrxnId} = $1;
      last;
    }
  }

  $rxn{$rxnId} = $obj;
}
close($fhr) || die "Error reading $indir/reactions.dat";
print STDERR "Read " . scalar(keys %rxn) . " reactions\n";

# EC numbers may be associated with the top-level reaction but not the subreactions.
# So, record this assocation so that later on the EC number can be "pushed down" to the subreactions.
my $nECDown = 0;
while (my ($rxnId, $rxn) = each %rxn) {
  foreach my $subId (@{ $rxn->{subreactions} }) {
    die "Invalid sub-reaction id $subId" unless exists $rxn{$subId};
    my $sub = $rxn{$subId};
    if (@{ $rxn->{ecnum} } > 0) {
      push @{ $sub->{ecnum} }, @{ $rxn->{ecnum} };
      $nECDown++;
    }
  }
}
print STDERR "Pushed $nECDown EC numbers from top-level reactions to sub-reactions\n";

open(OUT, ">", "$out/MetacycPathway.tab")
  || die "Cannot write to $out/MetacycPathway.tab";
foreach my $pathwayId (sort keys %path) {
  my $path = $path{$pathwayId};
  print OUT join("\t", $pathwayId, $path->{"pathwayName"})."\n";
}
close(OUT) || die "Error writing to $out/MetacycPathway.tab";

open(OUT, ">", "$out/MetacycPathwayReaction.tab")
  || die "Cannot write to $out/MetacycPathwayReaction.tab";
foreach my $pathwayId (sort keys %path) {
  my $path = $path{$pathwayId};
  # A super-pathway can include a reaction multiple times -- just ignore repeats
  my %seen = (); # rxnId => 1
  foreach my $rxn (@{ $path->{"reactions"} }) {
    my $rxnId = $rxn->{rxnId};
    next if exists $path{$rxnId};
    next if exists $seen{$rxnId};
    $seen{$rxnId} = 1;
    die "Pathway $pathwayId includes non-reaction non-pathway $rxnId"
      unless exists $rxn{$rxnId};
    print OUT join("\t", $pathwayId, $rxn->{"rxnId"}, $rxn->{"direction"},
                   exists $path->{hypo_reactions}{ $rxn->{rxnId} } ? 1 : 0)."\n";
  }
}
close(OUT) || die "Error writing to $out/MetacycPathwayReaction.tab";

open(OUT, ">", "$out/MetacycPathwayReactionPredecessor.tab")
  || die "Cannot write to $out/MetacycPathwayReactionPredecessor.tab";
foreach my $pathwayId (sort keys %path) {
  my $path = $path{$pathwayId};
  my $pred = $path->{predecessor};
  foreach my $rxnId (sort keys %{ $path->{predecessor} }) {
    next if exists $path{$rxnId};
    die "Pathway $pathwayId has predecessors for non-reaction non-pathway $rxnId"
      unless exists $rxn{$rxnId};
    # ensure that duplicate predecessors are not written (they are not meaningful)
    my %seen = (); # pre => 1
    foreach my $predId (sort @{ $path->{predecessor}{$rxnId} }) {
      next if exists $seen{$predId};
      $seen{$predId} = 1;
      print OUT join("\t", $pathwayId, $rxnId, $predId)."\n";
    }
  }
}
close(OUT) || die "Error writing to $out/MetacycPathwayReactionPredecessor.tab";

open(OUT, ">", "$out/MetacycPathwayPrimaryCompound.tab")
  || die "Cannot write to $out/MetacycPathwayPrimaryCompound.tab";
foreach my $pathwayId (sort keys %path) {
  my $path = $path{$pathwayId};
  my %seenRxn = (); # rxnId => 1 # ignore duplicates of reactions
  foreach my $rxn (@{ $path->{"reactions"} }) {
    my $rxnId = $rxn->{rxnId};
    next if exists $path{$rxnId};
    die "Pathway $pathwayId has non-reaction non-pathway component $rxnId"
      unless exists $rxn{$rxnId};
    next if exists $seenRxn{$rxnId};
    $seenRxn{$rxnId} = 1;
    # prevent duplicates, even though they may be deliberate in some biosynthetic pathways as a
    # hacked way to show the requirement for >1 molecule
    my %seen = (); # compound => side => 1
    foreach my $primary (@{ $rxn->{"primary"} }) {
      my ($compoundId, $side) = @$primary;
      next if exists $seen{$compoundId}{$side};
      print OUT join("\t", $pathwayId, $rxn->{"rxnId"}, $side, $compoundId)."\n";
      $seen{$compoundId}{$side} = 1;
    }
  }
}
close(OUT) || die "Error writing to $out/MetacycPathwayPrimaryCompound.tab";

open(OUT, ">", "$out/MetacycReaction.tab")
  || die "Cannot write to $out/MetacycReaction.tab";
foreach my $rxnId (sort keys %rxn) {
  my $rxn = $rxn{$rxnId};
  my $name = $rxn->{"rxnName"}
    || join(",", @{ $rxn->{"ecnum"} });
  print OUT join("\t", $rxnId, $name, $rxn->{isSpontaneous}, $rxn->{keggrxnId})."\n";
}
close(OUT) || die "Error writing to $out/MetacycReaction.tab";

open(OUT, ">", "$out/MetacycReactionCompound.tab")
  || die "Cannot write to $out/MetacycReactionCompound.tab";
foreach my $rxnId (sort keys %rxn) {
  my $rxn = $rxn{$rxnId};
  foreach my $cmp (@{ $rxn->{compounds} }) {
    print OUT join("\t", $rxnId,
                   $cmp->{"compoundId"}, $cmp->{"side"},
                   $cmp->{"coefficient"}, $cmp->{"compartment"}
                  )."\n";
  }
}

close(OUT) || die "Error writing to $out/MetacycReactionCompound.tab";

open(OUT, ">", "$out/MetacycReactionEC.tab")
  || die "Cannot write to $out/MetacycReactionEC.tab";
foreach my $rxnId (sort keys %rxn) {
  my $rxn = $rxn{$rxnId};
  my %seen = ();
  foreach my $ecnum (@{ $rxn->{"ecnum"} }) {
    next if exists $seen{$ecnum};
    $seen{$ecnum} = 1;
    print OUT join("\t", $rxnId, $ecnum)."\n" unless $ecnum =~ m/-/; # fully specified only
  }
}
close(OUT) || die "Error writing to $out/MetacycReactionEC.tab";

open(OUT, ">", "$out/MetacycCompound.tab")
  || die "Cannot write to $out/MetacycCompound.tab";
foreach my $compoundId (sort keys %compounds) {
  my $cmp = $compounds{$compoundId};
  my $name = $cmp->{"compoundName"};
  # Quote names for sqlite3's benefit
  if ($name =~ m/"/) {
    $name =~ s/"/""/g;
    $name = '"' . $name . '"';
  }
  print OUT join("\t", $compoundId,
                 $name, $cmp->{"keggLigand"}, $cmp->{"formula"},
                0 # not a class
                )."\n";
}
foreach my $id (sort keys %otherCompounds) {
  print OUT join("\t", $id, $classes{$id}{"name"}, "", "", 1) . "\n"
    if exists $classes{$id};
}
close(OUT) || die "Error writing to $out/MetacycCompound.tab";
