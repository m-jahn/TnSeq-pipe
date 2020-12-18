#!/usr/bin/perl -w
# Fetch information about a list of genes from the sqlite3 database
use strict;
use DBI;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()
use Utils; # for seed_desc

my $usage = <<END
fetchGenes.pl -db feba.db -genes geneId1,geneId,...geneId > out.tab
	or
fetchGenes.pl -db feba.db -in genesfile > out.tab

The output is tab-delimited with the fields
	orgId, locusId, sysName, original_description, SEED_description, KEGG_description

The input file may contain the fields orgId (or org), locusId, and sysName --
either of the latter two may suffice.

If specifying the genes on the command line, either orgId:locusId or sysName is acceptable.
END
;

my ($dbfile, $infile, $genespec);
die $usage if @ARGV == 0;
die $usage
    unless GetOptions('db=s' => \$dbfile,
                      'in=s' => \$infile,
                      'genes=s' => \$genespec)
    && @ARGV == 0;
die "Not all required arguments are present:\n$usage"
    unless defined $dbfile && (defined $infile || defined $genespec);
die "Cannot specify both -in and -genes" if defined $infile && defined $genespec;
die "Not a file: $dbfile\n" unless -e $dbfile;
die "Not a file: $infile\n" if defined $infile && ! -e $infile;

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 })
    || die $DBI::errstr;

my @genes = ();
if (defined $infile) {
    my @in = ReadTable($infile, []);
    die "No input rows\n" if scalar(@in) == 0;
    foreach my $row (@in) {
        my $orgId = $row->{orgId} || $row->{org};
        my $locusId = $row->{locusId};
        my $sysName = $row->{sysName};
        my @where = ();
        my @args = ();
        if ($orgId) {
            push @where, "orgId = ?";
            push @args, $orgId;
        }
        if ($sysName) {
            push @where, "sysName = ?";
            push @args, $sysName;
        } elsif ($locusId) {
            push @where, "locusId = ?";
            push @args, $locusId;
        } else {
            die "Must specify locusId or sysName";
        }
        die "No where clause" if @where == 0;
        my $where = join(" AND ", @where);
        my $query = "SELECT * from Gene WHERE $where";
        my $hits = $dbh->selectall_arrayref($query, { Slice => {} }, @args);
        die "No gene matches $where with arguments " . join(", ", @args) . "\n"
            if @$hits == 0;
        die "More than one gene matches $where with arguments " . join(", ", @args) . "\n"
            if @$hits > 1;
        push @genes, $hits->[0];
    }
} elsif (defined $genespec) {
    my @genespec = split /,/, $genespec;
    die "No genes specified\n" if @genespec == 0;
    foreach my $spec (@genespec) {
        my $where;
        my @args = ();
        if ($spec =~ /^(.*):(.*)$/) {
            my ($orgId,$name) = ($1,$2);
            $where = "orgId = ? AND (locusId = ? OR sysName = ?)";
            push @args, $orgId, $name, $name;
        } else {
            $where = "locusId = ? OR sysName = ?";
            push @args, $spec, $spec;
        }
        my $query = "SELECT * from Gene WHERE $where";
        my $hits = $dbh->selectall_arrayref($query, { Slice => {} }, @args);
        die "No gene matches $where with argument " . join(", ", @args) . "\n"
            if @$hits == 0;
        die "More than one gene matches $where with arguments " . join(", ", @args) . "\n"
            if @$hits > 1;
        push @genes, $hits->[0];
    }
} else {
    die;
}

# Remove duplicates
my %seen = ();
@genes = grep { my $id = $_->{orgId} . ":" . $_->{locusId};
                my $seen = exists $seen{$id};
                $seen{$id} = 1;
                ! $seen; } @genes;

print STDERR "Read " . scalar(@genes) . " genes\n";

print join("\t", qw{orgId locusId sysName specific_sick new_annotation original_description SEED_description KEGG_description})."\n";
foreach my $gene (@genes) {
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};
    my ($seed_desc, $seed_classes) = Utils::seed_desc($dbh, $orgId, $locusId);
    my $kegg_descs = $dbh->selectcol_arrayref("SELECT DISTINCT KgroupDesc.desc
                                               FROM BestHitKEGG JOIN KEGGMember USING (keggOrg,keggId)
                                               JOIN KgroupDesc USING (kgroup)
                                               WHERE orgId = ? AND locusId = ? AND desc <> ''",
                                              {}, $orgId, $locusId);
    my $spec = $dbh->selectall_arrayref(qq{ SELECT DISTINCT expGroup, condition FROM SpecOG
                                            WHERE orgId = ? AND locusId = ? AND minFit < 0 },
                                        { Slice => {} },
                                        $orgId, $locusId);
    my $reanno = $dbh->selectrow_hashref("SELECT * FROM Reannotation WHERE orgId = ? AND locusId = ?",
                                         {}, $orgId, $locusId);
    my @spec = map $_->{expGroup} . ": " . $_->{condition}, @$spec;
    print join("\t", $orgId, $locusId, $gene->{sysName},
               join("; ", @spec),
               $reanno->{new_annotation} || "",
               $gene->{desc},
               $seed_desc || "",
               join(": ", @$kegg_descs) || ""
        ) . "\n";
}

    
