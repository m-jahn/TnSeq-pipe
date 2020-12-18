#!/usr/bin/perl -w
# A helper script for db_setup.pl -- computes the SpecOG table
use strict;
use Getopt::Long;
use DBI;

my $db;
my $dir = ".";
my $noupdate;
my $usage = "Usage: db_setup_specOG.pl -db feba.db [ -dir $dir ] [ -noupdate ] [ -debug ]\n";
my $debug;

die $usage
    unless GetOptions('db=s' => \$db, 'dir=s' => \$dir, 'noupdate' => \$noupdate, 'debug' => \$debug)
    && defined $db
    && @ARGV == 0;

die "No such directory: $dir" unless -d $dir;
die "No such file: $db" unless -e $db;
my $dbh = DBI->connect("dbi:SQLite:dbname=$db", "", "", { RaiseError => 1 }) || die $DBI::errstr;

# First, metadata about all experiments that have specific phenotypes
my %specExp = (); # group => condition => orgId => expName => 1
my $expsData = $dbh->selectall_arrayref("SELECT DISTINCT orgId,expName,expGroup,condition_1
                                         FROM Experiment JOIN SpecificPhenotype USING (orgId,expName);",
                                         { Slice => {} });
foreach my $row (@$expsData) {
    $specExp{  $row->{expGroup} }{ $row->{condition_1} }{ $row->{orgId} }{ $row->{expName} } = 1;
}

# All specific phenotypes by gene
my %specGene = (); # group => condition => orgId => locusId => 1
my %orgGene = (); # orgId => locusId => 1
my $geneData = $dbh->selectall_arrayref("SELECT orgId,locusId,expName,expGroup,condition_1
                                         FROM Experiment JOIN SpecificPhenotype USING (orgId,expName);",
                                        { Slice => {} });
foreach my $row (@$geneData) {
    $orgGene{ $row->{orgId} }{ $row->{locusId} } = 1;
    $specGene{ $row->{expGroup} }{ $row->{condition_1} }{ $row->{orgId} }{ $row->{locusId} } = 1;
}

# All BBHs between all relevant genes
my %bbh = (); # org1 => locus1 => org2 => locus2
my $nPair = 0;
while (my ($orgId, $geneHash) = each %orgGene) {
    my @genes = keys(%$geneHash);
    if (@genes == 0) {
        print STDERR "Warning: no specific phenotypes in organism $orgId\n";
        next;
    }
    my $geneSpec = join(",", map {"'" . $_ . "'"} @genes);

    my $hits = $dbh->selectall_arrayref("SELECT locusId1,orgId2,locusId2 FROM Ortholog WHERE orgId1 = ?
                                        AND locusId1 IN ($geneSpec)",
                                       {}, $orgId);
    foreach my $row (@$hits) {
        my ($locus1,$org2,$locus2) = @$row;
        if (exists $orgGene{$org2}{$locus2}) {
            $bbh{$orgId}{$locus1}{$org2} = $locus2;
            $nPair++;
        }
    }
}

my $outfile = "$dir/db.SpecOG";
open(SPECOG, ">", $outfile) || die "Cannot write to $outfile";

# For each group and condition, cluster into OGs, fetch fitness values, and write out the results
my $nOG = 1;
while (my ($expGroup, $condHash) = each %specGene) {
    while (my ($condition, $orgHash) = each %$condHash) {
        # First, make a list of orgId/locusId pairs and sort by how many orthologs they have
        my @pairs = (); # orgId, locusId, nBBH
        while (my ($orgId, $geneHash) = each %$orgHash) {
            foreach my $locusId (keys %$geneHash) {
                push @pairs, [ $orgId, $locusId, scalar(keys %{ $bbh{$orgId}{$locusId} }) ];
            }
        }
        @pairs = sort { $b->[2] <=> $a->[2] } @pairs;
        # Go through the genes, assigning them to an OG if possible
        my %og = (); # orgId => locusId => nOG
        my %ogSize = (); # ogId => #members
        foreach my $pair (@pairs) {
            my ($orgId,$locusId,undef) = @$pair;
            unless (exists $og{$orgId}{$locusId}) {
                $og{$orgId}{$locusId} = $nOG;
                $ogSize{$nOG} = 1;
                $nOG++;
            }
            # and assign each BBH *that is relevant* to this OG
            my $bbhs = $bbh{$orgId}{$locusId};
            if (defined $bbhs) {
                while (my ($org2,$locus2) = each %{ $bbh{$orgId}{$locusId}}) {
                    next unless exists $orgHash->{$org2}{$locus2};
                    unless (exists $og{$org2}{$locus2}) {
                        $og{$org2}{$locus2} = $og{$orgId}{$locusId};
                        $ogSize{ $og{$orgId}{$locusId} }++;
                    }
                }
            }
        }
        # fetch the fitness values
        my %fit = (); # orgId => locusId => [ locusId, minfit, maxfit, minT, maxT ]
        while (my ($orgId, $geneHash) = each %$orgHash) {
            my @genes = keys %$geneHash;
            next unless @genes > 0; # empty hash may be created by exists check above
            my $geneSpec = join(",", map {"'".$_."'"} @genes);
            my @expNames = keys %{ $specExp{$expGroup}{$condition}{$orgId} };
            die unless @expNames > 0;
            my $expSpec = join(",", map {"'".$_."'"} @expNames);
            print STDERR "Group $expGroup condition $condition org $orgId geneSpec $geneSpec expSpec $expSpec\n"
              if defined $debug;
            my $rows = $dbh->selectall_arrayref(
                "SELECT locusId, min(fit), max(fit), min(t), max(t)
                 FROM GeneFitness
                 WHERE orgId = ? AND locusId IN ($geneSpec) AND expName IN ($expSpec)
                 GROUP BY locusId",
                {}, $orgId);
            foreach my $row (@$rows) {
                $fit{$orgId}{$row->[0]} = $row;
            }
        }
        # and write out the ogs
        while (my ($orgId, $geneHash) = each %$orgHash) {
            foreach my $locusId (keys %$geneHash) {
                die if !defined $og{$orgId}{$locusId};
                die "Locus $locusId has specific phenotype but no fitness data?"
                  unless defined $fit{$orgId}{$locusId};
                my (undef, $minfit, $maxfit, $minT, $maxT) = @{ $fit{$orgId}{$locusId} };
                die "No data for $orgId $locusId $expGroup $condition" unless defined $minfit;
                print SPECOG join("\t", $og{$orgId}{$locusId}, $expGroup, $condition,
                                  $orgId, $locusId,
                                  $minfit,$maxfit,$minT,$maxT,
                                  $ogSize{$og{$orgId}{$locusId}} )."\n";
            }
        }
    }
}
close(SPECOG) || die "Error writing to $outfile";

if (defined $noupdate) {
    print STDERR "Wrote to $outfile\n";
} else {
    print STDERR "Wrote temporary file: $outfile\n";
    $dbh->disconnect();
    open(SQLITE, "|-", "sqlite3", $db) || die "Cannot run sqlite3 on $db";
    print SQLITE <<END
.mode tabs
DELETE from SpecOG;
.import $outfile SpecOG
END
    ;
    close(SQLITE) || die "Error running sqlite3 on $db";
    print STDERR "Database updated, removing temporary file\n";
    unlink($outfile);
}
