#!/usr/bin/perl -w
# Create the metacyc_ids.tab and metacyc.faa files
# Tested with MetaCyc 22.6

# Convert the metacyc uniprot-seq-ids.dat into a simple tab-delimited format of
# rxn_id	ec_number	sequence_id
#
# Tested with MetaCyc 19.5 -- note this uses only what is in the metacyc specific database
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use ParsePTools;
use FEBA_Utils qw{ReadFastaEntry};

my $prefix = "metacyc";
my $prerapsearch = "$Bin/prerapsearch";

my $usage = <<END
Usage: ParseMetaCycSeqDat.pl -in metacyc_data_dir

Writes to metacyc_ids.tab (tab-delimited with a header line)
and metacyc.faa (in fasta format). Also formats the fasta database
using prerapsearch. Once this is run, RunRapSearch.pl can be
used to compute MetaCyc best hits, and db_setup_metacyc.pl
can be used to parse those hits (with the help of
metacyc_ids.tab). This script does not parse the
MetaCyc pathways; for that, see ParseMetaCycPathways.pl

To link protein sequences to reaction ids, uses CATALYZES entries for
proteins or for complexes, and then uses the REACTION entry for the
enzymatic reaction (which includes transport reactions). It
also gets an EC number from EC-NUMBER field of reactions.

Limitations -- multi-level enzyme complexes (A is a COMPONENT-OF B
which is a COMPONENT-OF C which has CATALYZES) are not handled.

Optional arguments:
  -out $prefix -- prefix of output files
  -prerapsearch $prerapsearch -- the prerapsearch executable
END
;

my $indir;
die $usage
  unless GetOptions('indir=s' => \$indir,
                    'prefix=s' => \$prefix,
                    'prerapsearch=s' => \$prerapsearch)
  && @ARGV == 0
  && defined $indir;
die "No such directory: $indir\n" unless -d $indir;
die "No such executable: $prerapsearch\n" unless -x $prerapsearch;

my @files = ("protseq.fsa", "proteins.dat", "enzrxns.dat", "reactions.dat");
foreach my $file (@files) {
  die "No such file: ${indir}/$file\n" unless -e "${indir}/$file";
}

# Save the sequence descriptions
my %desc = (); # monomer identifier to description
open(my $fhFaa, "<", "${indir}/protseq.fsa")
  || die "Cannot read ${indir}/protseq.fsa";
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fhFaa, $state)) {
  die "Empty sequence after $header" if $seq eq "";
  $header =~ m/^gnl[|]META[|](\S+) (.*)$/
    || die "Expected sequence named gnl|META|id but it is not found in this header: $header\n";
  my ($id, $desc) = ($1, $2);
  die "Duplicate entry for identifier $id in ${indir}/protseq.fsa\n"
    if exists $desc{$id};
  $desc =~ s/\s/ /g; # ensure no tabs or other unusual whitespace
  $desc{$id} = $desc;
}
close($fhFaa) || die "Error reading ${indir}/protseq.fsa";
print STDERR "Read " . scalar(keys %desc) . " sequences\n";

# Parse the proteins. Save CATALYZES and COMPONENT-OF relationships
my %catalyzes = (); # protein unique_id => list of enzymatic reactions
my %componentof = (); # protein unique_id => list of complexes' unique_ids
my %protId = (); # valid protein ids

open (my $fhProt, "<", "${indir}/proteins.dat")
  || die "Cannot read ${indir}/proteins.dat";
while (my $prot = ParsePTools($fhProt)) {
  my $id = $prot->{"UNIQUE-ID"}[0]{"value"};
  die unless $id;
  $protId{$id} = 1;
  if (exists $prot->{CATALYZES}) {
    my @catalyzes = map $_->{value}, @{ $prot->{CATALYZES} };
    $catalyzes{$id} = \@catalyzes;
  }
  if (exists $prot->{"COMPONENT-OF"}) {
    my @complexes = map $_->{value}, @{ $prot->{"COMPONENT-OF"} };
    $componentof{$id} = \@complexes;
  }
}
close($fhProt) || die "Error reading ${indir}/proteins.dat";
print STDERR "Read data for " . scalar(keys %protId) . " proteins and complexes\n";

# Verify that all proteins with sequence are known
my $nDirectLink = 0; # sequences that are linked to 1 or more enzrxns directly
foreach my $protId (keys %desc) {
  die "Monomer identifier $protId in ${indir}/protseq.fsa is not in proteins.dat\n"
    unless exists $protId{$protId};
  $nDirectLink++ if exists $catalyzes{$protId};
}

# Propagate enzrxn links to complex members
while (my ($protId, $complexes) = each %componentof) {
  foreach my $complexId (@$complexes) {
    push @{ $catalyzes{$protId} }, @{ $catalyzes{$complexId} }
      if exists $catalyzes{$complexId};
  }
}
my $nTotalLink = grep { exists $catalyzes{$_} } (keys %desc);
my $nComponentLink  = $nTotalLink - $nDirectLink;
my $nNoLink = scalar(keys %desc) - $nTotalLink;
print STDERR "Sequences linked to enzymatic reactions: $nDirectLink directly and $nComponentLink via complexes; $nNoLink not linked\n";

# Parse the enzymatic reactions. Save the link to reactions.
my %enzrxnToRxn = (); # enzrxn unique id to list of reaction ids
open(my $fhEnz, "<", "${indir}/enzrxns.dat")
  || die "Cannot read ${indir}/enzrxns.dat";
my $nEnzRxn = 0;
while (my $enzrxn = ParsePTools($fhEnz)) {
  my $id = $enzrxn->{"UNIQUE-ID"}[0]{"value"};
  die unless $id;
  $nEnzRxn++;
  if (exists $enzrxn->{REACTION}) {
    my @rxnIds = map $_->{value}, @{ $enzrxn->{REACTION} };
    $enzrxnToRxn{$id} = \@rxnIds;
  }
}
close($fhEnz) || die "Error reading ${indir}/enzrxns.dat";
print STDERR "Read data for $nEnzRxn enzymatic reactions with " . scalar(keys %enzrxnToRxn) . " linked to reactions\n";

# Parse the reactions. Save link to EC numbers.
# Also save that the reaction is valid (as an empty EC list).
my %rxnEc = (); # reaction to list of EC numbers (or empty)
open(my $fhRxn, "<", "${indir}/reactions.dat")
  || die "Cannot read ${indir}/reactions.dat";
while (my $rxn = ParsePTools($fhRxn)) {
  my $rxnId = $rxn->{"UNIQUE-ID"}[0]{"value"};
  die unless $rxnId;
  my @ecs = ();
  @ecs = map $_->{value}, @{ $rxn->{"EC-NUMBER"}}
    if exists $rxn->{"EC-NUMBER"};
  @ecs = map { s/^EC-//; $_ } @ecs; # MetaCyc entries usually have this prefix
  $rxnEc{$rxnId} = \@ecs;
}
close($fhRxn) || die "Error reading ${indir}/reactions.dat";
print STDERR "Read data for " . scalar(keys %rxnEc) . " reactions\n";

my $outtab = "${prefix}_ids.tab";
open(my $fhTab, ">", $outtab) || die "Cannot write to $outtab\n";
# May report multiple lines for one hit
print $fhTab join("\t", qw{protId desc rxnId EC})."\n";
foreach my $protId (sort keys %desc) {
  my $desc = $desc{$protId};
  my $enzrxns = [];
  $enzrxns = $catalyzes{$protId} if exists $catalyzes{$protId};
  my @rxnIds = ();
  foreach my $enzrxnId (@$enzrxns) {
    push @rxnIds, @{ $enzrxnToRxn{$enzrxnId} } if exists $enzrxnToRxn{$enzrxnId};
  }
  # eliminate any duplicates
  my %sofar = ();
  @rxnIds = grep { my $keep = !exists $sofar{$_}; $sofar{$_} = 1; $keep } @rxnIds;
  if (@rxnIds > 0) {
    foreach my $rxnId (@rxnIds) {
      my $ec = exists $rxnEc{$rxnId} ? $rxnEc{$rxnId} : [];
      if (@$ec > 0) {
        foreach my $ec (@$ec) {
          print $fhTab join("\t", $protId, $desc, $rxnId, $ec)."\n";
        }
      } else {
        print $fhTab join("\t", $protId, $desc, $rxnId, "")."\n";
      }
    }
  } else {
    print $fhTab join("\t", $protId, $desc, "", "")."\n";
  }
}
close($fhTab) || die "Error writing to $outtab";
print STDERR "Wrote $outtab\n";

my $outseq = "$prefix.faa";
system("cp", "${indir}/protseq.fsa", $outseq) == 0
  || die "Error copying to $outseq\n";
system($prerapsearch, "-d", $outseq, "-n", "$prefix.rap", "-f", "T") == 0
|| die "Error running prerapsearch\n";
print STDERR "Created $prefix.faa and $prefix.rap\n";
