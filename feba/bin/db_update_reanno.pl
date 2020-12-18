#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()
use DBI;

my $usage = <<END
Usage: db_setup.pl -reanno reannotation_file -db db_file_name

    reannotation_file should be tab-delimited and should include the
    fields orgId, locusId, new_annotation, and comment.

Optional arguments:
  -nometacyc -- skip recomputing metacyc pathway coverage
  -debug -- report progress
END
    ;

{
    my ($reannofile, $db, $nometacyc, $debug);
    die $usage unless GetOptions('reanno=s' => \$reannofile,
                                 'db=s' => \$db,
                                 'nometacyc' => \$nometacyc,
                                 'debug' => \$debug)
        && @ARGV == 0
        && defined $reannofile
        && defined $db;
    my @annos = ReadTable($reannofile, ['orgId', 'locusId', 'new_annotation', 'comment']);
    print STDERR "Read " . scalar(@annos) . " rows from $reannofile\n";
    my $tmpDb = "$db.$$.tmp";
    print STDERR "Copying to $tmpDb (slow)\n";
    system("cp", $db, $tmpDb) == 0 || die "Cannot cp $db $tmpDb -- $!";

    my $dbh = DBI->connect("dbi:SQLite:dbname=$tmpDb", "", "", { RaiseError => 1 }) || die $DBI::errstr;
    print STDERR "Modifying $tmpDb\n";

    my $clear_statement = $dbh->prepare(qq{ DELETE FROM Reannotation; });
    $clear_statement->execute() || die "Cannot clear the Reannotation table in $tmpDb";
    my $clear_statement2 = $dbh->prepare(qq{ DELETE FROM ReannotationEC; });
    $clear_statement2->execute() || die "Cannot clear the ReannotationEC table in $tmpDb";

    my $insert_statement = $dbh->prepare(qq{ INSERT INTO Reannotation VALUES (?,?,?,?); });
    my $insert_statement2 = $dbh->prepare(qq{ INSERT INTO ReannotationEC VALUES (?,?,?); });

    my $orgs = $dbh->selectcol_arrayref("SELECT orgId FROM Organism");
    my %orgs = map { $_ => 1 } @$orgs;
    my $nAdd = 0;
    my %skipOrgs = ();
    my @ecrows = (); # orgId, locusId, ec
    foreach my $row (@annos) {
      if (! exists $orgs{ $row->{orgId} }) {
        $skipOrgs{ $row->{orgId} } = 1;
      } else {
        my ($n) = $dbh->selectrow_array("SELECT COUNT(*) FROM Gene WHERE orgId = ? AND locusId = ?",
                                        {}, $row->{orgId}, $row->{locusId});
        if ($n == 0) {
          print STDERR "Warning: reannotation of non-existent gene $row->{locusId} in $row->{orgId}\n";
        } else {
          $insert_statement->execute($row->{orgId}, $row->{locusId}, $row->{new_annotation}, $row->{comment})
            || die "Failed to insert into Reannotation";
          print "Inserted $row->{orgId} $row->{locusId}\n" if $debug;
          $nAdd++;
          # And try to parse EC number(s) from the annotation
          # This can be [EC 1.1.1.1] or [EC 1.1.1.1 1.1.1.2] or [EC 1.1.1.1; 1.1.1.1.2]
          # or similarly with commas or parens, or it can be
          # (EC 1.1.1.1; EC 1.1.1.1)
          my @words = split / /, $row->{new_annotation};
          my %ecSeen = ();
          foreach my $word (@words) {
            # Removing leading paren or bracket
            $word =~ s/^[(\[]//;
            # Remove leading EC: (the space case does not arise)
            $word =~ s/^EC://i;
            $word =~ s/[,;)\]]+$//;
            # allow EC numbers like 3.2.1.B4, sucrose-6-phosphate hydrolase
            if ($word =~ m/^\d+[.]\d+[.]\d+[.][A-Z]?[0-9-]+$/) {
              my $ec = $word;
              if (exists $ecSeen{$ec}) {
                print STDERR "Warning: duplicate ec assignment in $row->{orgId} $row->{locusId} $row->{new_annotation}\n";
              } else {
                $insert_statement2->execute($row->{orgId}, $row->{locusId}, $ec)
                  || die "Failed to insert into ReannotationEC";
                $ecSeen{$ec} = 1;
              }
            } elsif ($word =~ m/\d+[.]\d+[.]\d+/) {
              print STDERR "Warning: unparseable ec number in $row->{orgId} $row->{locusId} $row->{new_annotation}\n";
            }
          }
        }
      }
    }
    $dbh->disconnect();
    print STDERR "Cleared Reannotation table and added $nAdd rows\n";
    print STDERR "Ignored reannotations for unknown organisms: " . join(" ", sort keys %skipOrgs) . "\n"
      if keys(%skipOrgs) > 0;
    unless (defined $nometacyc) {
      my @covercmd = ("$Bin/db_update_metacyc_coverage.pl", "-db", $tmpDb);
      system(@covercmd) == 0 || die "Error running " . join(" ", @covercmd) . "\n: $!";
    }

    print STDERR "Moving $tmpDb to $db\n" if $debug;
    rename($tmpDb, $db) || die "Cannot rename $tmpDb to $db\n";
}
