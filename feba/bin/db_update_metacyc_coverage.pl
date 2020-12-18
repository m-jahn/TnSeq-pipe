#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Utils; # for EcGenesToAll

my $usage = <<END
Usage: db_update_metacyc_coverage.pl -db db_file_name

   Empties and updates the MetacycPathwayCoverage table in the database
END
;

{
  my ($db);
  die $usage unless GetOptions('db=s' => \$db)
    && @ARGV == 0
    && defined $db;
  my $dbh = DBI->connect("dbi:SQLite:dbname=$db", "", "", { RaiseError => 1 }) || die $DBI::errstr;

  # For each reaction, find potential genes
  # using BestHitMetacyc, EC numbers, or SEED annotations (via links to KEGG reactions)
  my %rxnFound = (); # rxnId => orgId => 1 if found

  # By best hit
  my $rxn2locus = $dbh->selectall_arrayref("SELECT rxnId, orgId, locusId FROM BestHitMetacyc");
  foreach my $row (@$rxn2locus) {
    my ($rxnId,$orgId,$locusId) = @$row;
    $rxnFound{$rxnId}{$orgId} = 1;
  }

  # by EC numbers
  my $rxn2ec = $dbh->selectall_arrayref("SELECT rxnId, ecnum FROM MetacycReactionEC");
  my $orgs = $dbh->selectcol_arrayref("SELECT orgId FROM Organism");

  foreach my $orgId (@$orgs) {
    my $ecGenes = Utils::EcToGenesAll($dbh, $orgId); # ec => locusId => 1

    foreach my $row (@$rxn2ec) {
      my ($rxnId,$ec) = @$row;
      $rxnFound{$rxnId}{$orgId} = 1 if exists $ecGenes->{$ec};
    }
  }

  # by SEED
  foreach my $orgId (@$orgs) {
    my $rxnFoundSEED = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT rxnId FROM SEEDAnnotation
                                                    JOIN SEEDAnnotationToRoles USING (seed_desc)
                                                    JOIN SEEDRoleReaction USING (seedrole)
                                                    JOIN SEEDReaction USING (seedrxnId)
                                                    JOIN MetacycReaction USING (keggrxnId)
                                                    WHERE orgId = ? AND keggrxnId <> "" }, {}, $orgId);
    foreach my $rxnId (@$rxnFoundSEED) {
      $rxnFound{$rxnId}{$orgId} = 1;
    }
  }

  print STDERR "Found at least one candidate gene for " . scalar(keys %rxnFound) . " reactions\n";

  my $rxnSpont = $dbh->selectcol_arrayref("SELECT rxnId FROM MetacycReaction WHERE isSpontaneous = 1");
  # rxnId => 1 if spontaenous
  my %rxnSpont = map { $_ => 1 } @$rxnSpont;

  my $pathrxn = $dbh->selectall_arrayref("SELECT pathwayId, rxnId FROM MetacycPathwayReaction");
  my %pathRxns = ();
  foreach my $row (@$pathrxn) {
    my ($pathId,$rxnId) = @$row;
    push @{ $pathRxns{$pathId} }, $rxnId;
  }

  # Write a temporary file with MetaCycPathwayCoverage (orgId, pathwayId, nSteps, nFound)
  my $tmpdir = $ENV{TMPDIR} || "/tmp";
  my $tmpfile = "$tmpdir/db.MetacycPathwayCoverage.$$";
  my $nOut = 0;
  open(OUT, ">", $tmpfile) || die "Error writing to $tmpfile";
  foreach my $orgId (sort @$orgs) {
    foreach my $pathId (sort keys %pathRxns) {
      my @rxns = @{ $pathRxns{$pathId} };
      my $nSteps = scalar(@rxns);
      my $nWithGene = 0;
      my $nFound = 0;
      foreach my $rxnId (@rxns) {
        $nWithGene++ if $rxnFound{$rxnId}{$orgId};
        $nFound++ if $rxnFound{$rxnId}{$orgId} || exists $rxnSpont{$rxnId};
      }
      if ($nWithGene > 0) {
        print OUT join("\t", $orgId, $pathId, $nSteps, $nFound)."\n";
        $nOut++;
      }
    }
  }
  close(OUT) || die "Error writing to $tmpfile\n";
  print STDERR "Wrote temporary file $tmpfile\n";

  $dbh->disconnect();

  open(SQLITE, "|-", "sqlite3", $db) || die "Cannot run sqlite3 on $db";
  print SQLITE <<END
.mode tabs
DELETE FROM MetacycPathwayCoverage;
.import $tmpfile MetacycPathwayCoverage
END
;
  close(SQLITE) || die "Error running sqlite3 on $db";
  unlink($tmpfile);
  print STDERR "Updated MetacycPathwayCoverage in $db with $nOut rows\n";
}
