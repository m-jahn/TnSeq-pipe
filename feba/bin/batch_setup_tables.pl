#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadFastaDesc()

my $usage = <<END
Usage: batch_setup_tables.pl -db feba.db -aaseqs aaseqs -dir jobdir

The jobdir is assumed to contain two files:
	query.faa (queries)
	hits (usearch output)

The aaseqs has the proteins that queries were compared to.

The feba.db is the main database for the CGI scripts
(which is produced by db_setup.pl)

Creates the job.db database in that directory
END
    ;

my $minCoverage = 0.75;
my $minIdentity = 25.0;
my $minCofit = 0.8;

my $create_table_sql = <<END

CREATE TABLE Query(
    queryId TEXT NOT NULL,
    desc TEXT NOT NULL,
    isHypo TINYINT NOT NULL,
    aaseq TEXT NOT NULL,
    hasSpecific TINYINT NOT NULL,
    hasSpecificStrong TINYINT NOT NULL,
    hasSpecificCons TINYINT NOT NULL,
    hasCofit TINYINT NOT NULL,  /* always set if hasCofitCons is set, even if is below $minCofit */
    hasCofitCons TINYINT NOT NULL,
    PRIMARY KEY(queryId)
);

CREATE TABLE BestHit(
    queryId TEXT NOT NULL,
    orgId TEXT NOT NULL,
    locusId TEXT NOT NULL,
    identity REAL NOT NULL,
    coverage REAL NOT NULL,
    bits REAL NOT NULL,
    isBBH TINYINT NOT NULL,
    PRIMARY KEY (queryId,orgId,locusId)
);

END
    ;


{
    my ($dbfile, $aaseqFile, $dir);
    (GetOptions('db=s' => \$dbfile,
                'aaseqs=s' => \$aaseqFile,
                'dir=s' => \$dir)
     && defined $dbfile && defined $dir && defined $aaseqFile && @ARGV==0)
        || die $usage;

    die "No such file: $aaseqFile\n" unless -e $aaseqFile;
    die "No query.faa file in $dir\n" unless -e "$dir/query.faa";
    die "No hits in $dir\n" unless -e "$dir/hits";

    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","",{ RaiseError => 1 })
        || die $DBI::errstr;
    print STDERR "Analyzing hits in $dir/ against $dbfile\n";
    my %qinfo = ReadFastaDesc("$dir/query.faa");
    die $qinfo{'error'} if exists $qinfo{'error'};
    my $qseq = $qinfo{'seq'};
    my $qdesc = $qinfo{'desc'};
    my $qlen = $qinfo{'len'};

    my %sinfo = ReadFastaDesc($aaseqFile);
    die $sinfo{'error'} if exists $sinfo{'error'};
    my $slen = $sinfo{'len'};

    print STDERR "Parsed " . scalar(keys %$qdesc) . " query sequences and " . scalar(keys %$slen) . " subject sequences\n";

    my %besthit = (); # query => org => [ locusId, identity, coverage, bits ]
    my ($LOCUSID,$IDENTITY,$COVERAGE,$BITS) = 0..3;
    my %bestrev = (); # org => locusId => [ queryId, bits ]
    my $nHits = 0;
    open(HITS, "<", "$dir/hits") || die "Cannot read $dir/hits";
    while(<HITS>) {
        chomp;
        my ($queryId, $subjectId, $identity, $alnlen, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = split /\t/, $_, -1;
        die "Cannot parse ublast line:\n$_\n"
            unless defined $bits
            && $qBeg =~ m/^\d+$/ && $qEnd =~ m/^\d+$/
            && $sBeg =~ m/^\d+$/ && $sEnd =~ m/^\d+$/;
        die "Invalid queryId $queryId" unless exists $qlen->{$queryId};
        die "Invalid subjectId $subjectId" unless exists $slen->{$subjectId};
        $nHits++;
        my ($orgId,$locusId) = split /:/, $subjectId;
        die "Invalid subject $subjectId" unless defined $orgId && defined $locusId && $locusId ne "";

        unless (exists $bestrev{$orgId}{$locusId} && $bestrev{$orgId}{$locusId}[1] >= $bits) {
            $bestrev{$orgId}{$locusId} = [ $queryId, $bits ];
        }

        next if exists $besthit{$queryId}
        && exists $besthit{$queryId}{$orgId} 
        && $besthit{$queryId}{$orgId}[$BITS] >= $bits;

        my $cov1 = ($qEnd-$qBeg+1) / $qlen->{$queryId};
        my $cov2 = ($sEnd-$sBeg+1) / $slen->{$subjectId};
        my $cov = $cov1 < $cov2 ? $cov1 : $cov2;
        $besthit{$queryId}{$orgId} = [ $locusId, $identity, $cov, $bits ];
    }
    close(HITS) || die "Error reading $dir/hits";
    print STDERR "Parsed $nHits hits\n";

    # Load all specific phenotypes
    my $specs = $dbh->selectall_arrayref("SELECT * FROM SpecOG", { Slice => {} });
    my %spec = (); # orgId => locusId => list of [ expGroup, condition, minFit, maxFit, minT, maxT, nInOG ]
    foreach my $row (@$specs) {
        push @{ $spec{ $row->{orgId} }{ $row->{locusId} } },
        [ $row->{expGroup}, $row->{condition}, $row->{minFit}, $row->{maxFit}, $row->{minT}, $row->{maxT}, $row->{nInOG} ]; 
    }
    print STDERR "Loaded specific phenotypes\n";

    # Load all cofitness above $minCofit
    my $cofit = $dbh->selectall_arrayref("SELECT orgId,locusId FROM Cofit WHERE rank=1 AND cofit >= ?",
                                         {}, $minCofit);
    my %hascofit = (); # orgId => locusId => 1 if it has cofitness above threshold
    foreach my $row (@$cofit) {
        my ($orgId,$locusId) = @$row;
        $hascofit{$orgId}{$locusId} = 1;
    }
    print STDERR "Loaded cofitness\n";

    my $conscofit = $dbh->selectall_arrayref("SELECT DISTINCT orgId,locusId FROM ConservedCofit",
                                             {});
    my %hascofitCons = (); # org => locusId => 1 if conserved
    # Note -- the database uses a lower threshold for conserved cofitness, so
    # hasconscofit{orgId}{locusId} may be true even if hascofit{orgId}{locusId} is not
    foreach my $row (@$conscofit) {
        my ($orgId,$locusId) = @$row;
        $hascofitCons{$orgId}{$locusId} = 1;
    }

    # Report all the good hits, and report if they are bbhs
    # And record functional links
    my %qspec = (); # queryId => 1 if has specific
    my %qspecStrong = ();
    my %qspecCons = ();
    my %qcofit = (); # queryId => 1 if has cofit
    my %qcofitCons = ();
    my $bhFile = "$dir/db.BestHit";
    my $nBestHit = 0;
    my $nBBH = 0;
    open(OUT, ">", $bhFile) || die "Cannot write to $bhFile";
    foreach my $queryId (sort keys %besthit) {
        while (my ($orgId,$row) = each %{ $besthit{$queryId} }) {
            my ($locusId,$identity,$coverage,$bits) = @$row;
            next unless $identity >= $minIdentity && $coverage >= $minCoverage;
            $nBestHit++;
            my $bbh = $bestrev{$orgId}{$locusId}[0] eq $queryId ? 1 : 0;
            $nBBH++ if $bbh;
            print OUT join("\t", $queryId, $orgId, $locusId, $identity, $coverage, $bits, $bbh)."\n";
            if ($bbh) {
                if (exists $spec{$orgId}{$locusId}) {
                    $qspec{$queryId} = 1;
                    foreach my $row (@{ $spec{$orgId}{$locusId} }) {
                        my ($group,$cond,$minfit,$maxfit,$minT,$maxT,$nInOG) = @$row;
                        $qspecStrong{$queryId} = 1 if $minfit <= -2;
                        $qspecCons{$queryId} = 1 if $nInOG > 1;
                    }
                }
                $qcofit{$queryId} = 1 if exists $hascofit{$orgId}{$locusId};
                if (exists $hascofitCons{$orgId}{$locusId}) {
                    $qcofit{$queryId} = 1; # always highlight conserved cofitness
                    $qcofitCons{$queryId} = 1;
                }
            }
        }
    }
    close(OUT) || die "Error writing to $bhFile";
    print STDERR "Wrote $nBestHit best hits ($nBBH bidirectional) to $bhFile\n";

    my $qFile = "$dir/db.Query";
    open(OUT, ">", $qFile) || die "Cannot write to $qFile";
    foreach my $queryId (sort keys %$qseq) {
        my $desc = $qdesc->{$queryId};
        my $hypo = $desc eq ""
            || $desc =~ m/unknown function/i
            || $desc =~ m/uncharacterized/i
            || $desc =~ m/hypothetical/i
            || $desc =~ m/hypothtical/i
            || $desc =~ m/family/i
            || $desc =~ m/domain protein/i
            || $desc =~ m/related protein/i
            || $desc =~ m/transporter related/i;
        print OUT join("\t", $queryId, $desc,
                       $hypo ? 1 : 0,
                       $qseq->{$queryId},
                       $qspec{$queryId} ? 1 : 0,
                       $qspecStrong{$queryId} ? 1 : 0,
                       $qspecCons{$queryId} ? 1 : 0,
                       $qcofit{$queryId} ? 1 : 0,
                       $qcofitCons{$queryId} ? 1 : 0 )."\n";
    }
    close(OUT) || die "Error writing to $qFile";
    print STDERR "Wrote " . scalar(keys %$qseq) . " queries to $qFile\n";

    $dbh->disconnect();

    print STDERR "Loading tables\n";
    my $tmpdbfile = "$dir/tmp.$$.db";
    open(SQLITE, "|-", "sqlite3", "$tmpdbfile") || die "Cannot run sqlite3 on $tmpdbfile";
    print SQLITE ".bail on\n";
    print SQLITE ".mode tabs\n";
    print SQLITE "$create_table_sql\n";
    print SQLITE ".import $qFile Query\n";
    print SQLITE ".import $bhFile BestHit\n";
    close(SQLITE) || die "Error running sqlite3 commands\n";

    # and double check #rows in BestHit, the last table loaded
    my $dbh2 = DBI->connect("dbi:SQLite:dbname=$tmpdbfile", "", "", { RaiseError => 1 }) || die $DBI::errstr;
    my ($nBHActual) = $dbh2->selectrow_array("SELECT COUNT(*) FROM BestHit");
    die "Counting tables in BestHit failed" unless defined $nBHActual;
    die "Failed to load BestHit: expected $nBestHit but see $nBHActual instead"
        unless $nBHActual == $nBestHit;
    $dbh2->disconnect();

    rename($tmpdbfile, "$dir/job.db") || die "Cannot rename $tmpdbfile to $dir/job.db";
    print STDERR "Success\n";
}
