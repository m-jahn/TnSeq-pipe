#!/usr/bin/perl -w
# conserved_cofit.pl -- given an orthologs table and multiple fitness tables, identify
# pairs of genes & orthologs with conserved cofitness
use strict;
use Getopt::Long;

sub ReadOrthTable($);
sub ReadCofitFile($$$$$$);

{
    my $usage = qq{
Usage: conserved_cofit.pl -orth ortholog_table -out prefix [ -rank 10 -cor 0.6 ] [ -table ]
                          org1 cofitfile1 ... orgN cofitfileN

The ortholog table should either be a simple table with
tax1,locus1,tax2,locus2,ratio, or a matrix with locusId, taxonomyId,
and orthNNNN.locusId, where NNNN is the taxonomyId and orthologs in
its own genome are ignored. In either case, locusIds are of the form
organism:geneId.

The cofitness files (1 per organism) should be tab-delimited with no
org: prefix in the locusId or hitId. They should either have a header
line, or have fields in order orgId, locusId, hitId, rank, cofit.
(orgId is not needed if there is a header line.)

By default, all pairs are output to prefix.pairs, but the rank and R
options can be used to set thresholds for maximum rank and minimum
cofitness.

With the -table option, the file is suitable for loading into the database.
There is no header line and the pairs are output both ways. (By default,
only lines with locus1 < locus2 are shown.)

The order of fields in the output is: taxId locus1 locus2 rank cofit
otax orth1 orth2 orank ocofit.
};
    my $orthfile;
    my $out;
    my $maxrank;
    my $mincor;
    my $tableMode;
    (GetOptions('orth=s' => \$orthfile, 'out=s' => \$out,
                'table' => \$tableMode,
                'rank=i' => \$maxrank, 'cor=f' => \$mincor)
     && defined $orthfile && defined $out) || die $usage;
    die "No cofitness files" unless @ARGV > 0;
    die "Must have pairs of organism cofitness_file\n" unless scalar(@ARGV) % 2 == 0;
    my %cofitfiles = @ARGV; # organism => file

    # taxa: taxonomyId or nickname => locusId => 1
    # orth: locusId => taxonomyId => orthologId if it exists
    my ($taxa,$orth) = ReadOrthTable($orthfile);
    print STDERR "Read orthologs for " . scalar(keys %$orth) . " genes in ".scalar(keys %$taxa)." genomes\n";

    my %nTaxWithOrth = (); # for each taxon, the number of genes that have at least one ortholog
    foreach my $taxonomyId (sort keys %$taxa) {
	my $nWithOrth = 0;
	my $nGene = 0;
	foreach my $locusId (keys %{ $taxa->{$taxonomyId} }) {
	    $nGene++;
	    $nWithOrth++ if scalar(keys(%{ $orth->{$locusId} })) > 0;
	}
	print STDERR join("\t","nWithOrth",$taxonomyId,$nGene,$nWithOrth)."\n";
    }

    my %cofit = (); # taxon => locusId => hitId => [rank,cofitness]
    while (my ($taxonomyId, $cofitfile) = each %cofitfiles) {
	die "More than one cofitness file for $taxonomyId" if exists $cofit{$taxonomyId};
        # Only loads those with low rank and high correlation, if those options were used
	$cofit{$taxonomyId} = ReadCofitFile($cofitfile, $taxonomyId, $taxa, $orth, $maxrank, $mincor);
    }
    foreach my $taxonomyId (sort keys %$taxa) {
	print STDERR "No cofitness file for $taxonomyId\n" unless exists $cofit{$taxonomyId};
    }
    print STDERR "Read " . scalar(keys %cofit) . " cofitness tables\n";

    # And find pairs that have hits both ways in two organisms
    # Orthologs are not necessarily bidirectional or transitive, so do not handle that aspect of it
    my $nFound = 0;
    open(OUT,">","$out.pairs") || die "Cannot write to $out.pairs";
    print OUT join("\t",qw{taxId locus1 locus2 rank cofit otax orth1 orth2 orank ocofit})."\n"
        unless defined $tableMode;

    foreach my $taxonomyId (sort keys %$taxa) {
	# gene1 => gene2 => taxon => 1 if there is conserved cofitness in that organism
	# only direction gene1 < gene2 is stored
	my %pairs = ();
	my @othertax = grep { $_ ne $taxonomyId } (keys %$taxa);
	while (my ($locus1, $hithash) = each %{ $cofit{$taxonomyId} }) {
	    while (my ($locus2, $hitrow) = each %$hithash) {
		next unless defined $tableMode || ($locus1 cmp $locus2) < 0;
		my ($rank1,$cofit1) = @$hitrow;
                # Note -- if using minrank, this forces it to meet the rank threshold both ways
		next unless exists $cofit{$taxonomyId}{$locus2}{$locus1}; # must hit both ways
		my ($rank2,$cofit2) = @{ $cofit{$taxonomyId}{$locus2}{$locus1} };
		my $rank = $rank1 > $rank2 ? $rank1 : $rank2; # pessimistic estimate of importance
		die "cofitness for $locus1 $locus2 in taxon $taxonomyId is inconsistent: $cofit1 vs. $cofit2"
		    unless abs($cofit1 - $cofit2) <= 0.001;
		foreach my $othertax (@othertax) {
		    my $orth1 = $orth->{$locus1}{$othertax};
		    my $orth2 = $orth->{$locus2}{$othertax};
		    if (defined $orth1 && defined $orth2
			&& exists $cofit{$othertax}{$orth1}{$orth2} && exists $cofit{$othertax}{$orth2}{$orth1}) {
			$pairs{$locus1}{$locus2}{$othertax} = 1;
			$nFound++;
			my $orank1 = $cofit{$othertax}{$orth1}{$orth2}[0];
			my $orank2 = $cofit{$othertax}{$orth1}{$orth2}[0];
			my $orank = $orank1 > $orank2 ? $orank1 : $orank2;
			my $ocofit = $cofit{$othertax}{$orth1}{$orth2}[1];
			my (undef,$locus1Show) = split /:/, $locus1;
			my (undef,$locus2Show) = split /:/, $locus2;
			my (undef,$orth1Show) = split /:/, $orth1;
			my (undef,$orth2Show) = split /:/, $orth2;
			print OUT join("\t",$taxonomyId,$locus1Show,$locus2Show,$rank,$cofit1,
				       $othertax,$orth1Show,$orth2Show,$orank,$ocofit)."\n";
		    }
		}
	    }
	}
    }
    close(OUT) || die "Error writing to $out.pairs";
    print STDERR "Found $nFound cases of conserved cofitness (counting orthologs both ways)\n";
    print STDERR "Wrote $out.pairs\n";
}

# Returns references to taxa and orth hashes
sub ReadOrthTable($) {
    my ($orthfile) = @_;
    open(ORTH, "<", $orthfile) || die "Cannot read $orthfile";
    my %taxa = (); # taxonomyId => locusId => 1
    my %orth = ();

    my $header = <ORTH>;
    chomp $header;
    my @colnames = split /\t/, $header;
    my @expected = qw{tax1 locus1 tax2 locus2 ratio};
    if ($colnames[0] eq $expected[0] && scalar(@colnames) == scalar(@expected)) {
        # reading a .scores file with 5 fields
        foreach my $i (0..(scalar(@expected)-1)) {
            die "Column name in ortholog table is $colnames[$i] not $expected[$i]"
                unless $colnames[$i] eq $expected[$i];
        }
        while(my $line = <ORTH>) {
            chomp $line;
            my @F = split /\t/, $line, -1;
            die "Wrong number of fields in ortholog table:\n@F" unless scalar(@F) == scalar(@expected);
            my ($tax1,$locus1,$tax2,$locus2,undef) = @F;
            die "locus1 $locus1 has no :" unless $locus1 =~ m/:/;
            die "locus2 $locus2 has no :" unless $locus2 =~ m/:/;
            $taxa{$tax1}{$locus1} = 1;
            $taxa{$tax2}{$locus2} = 1;
            die "Duplicate entry for BBH of $locus1 in $tax2" if exists $orth{$locus1}{$tax2};
            $orth{$locus1}{$tax2} = $locus2;
        }
    } else {
        # reading the matrix style input file
        die "First column in ortholog table should be locus name" unless shift(@colnames) eq "locusId";
        die "Second column in ortholog table should be locus name" unless shift(@colnames) eq "taxonomyId";
        my @taxa = (); # in the order in the column names
        foreach my $name (@colnames) {
            $name =~ m/^orth(.*)[.][a-zA-z]+$/ || die "invalid ortholog field $name";
            push @taxa, $1;
        }
        # taxa: taxonomyId => locusId => 1
        %taxa = map { $_ => {} } @taxa;
        die "Duplicate taxon names" if scalar(keys %taxa) != scalar(@taxa);
        while(my $line = <ORTH>) {
            chomp $line;
            my @F = split /\t/, $line, -1;
            die "Wrong number of fields in ortholog table:\n@F" unless scalar(@F) == scalar(@taxa)+2;
            my $locusId = shift @F;
            die "Duplicate entry for locus $locusId" if exists $orth{$locusId};
            $orth{$locusId} = {}; # so we know it is a valid locusId
            my $taxonomyId = shift @F;
            $taxa{$taxonomyId}{$locusId} = 1;
            die "Invalid taxonomy $taxonomyId in ortholog table" unless exists $taxa{$taxonomyId};
            foreach my $orthtax (@taxa) {
                my $orthId = shift @F;
                die $line if !defined $orthId;
                $orth{$locusId}{$orthtax} = $orthId if $orthId ne "" && $orthtax ne $taxonomyId;
            }
        }
    }
    close(ORTH) || die "Error reading $orthfile";
    return(\%taxa, \%orth);
}

# inputs: filename, taxonomyId,
#	hash of taxa => taxonomyId => 1,
#	hash of locusId => taxonomyId => orthologId,
#	maxrank (or undef), mincor (or undef)
# returns a hash of locusId => hitId => [rank,cofitness]

sub ReadCofitFile($$$$$$) {
    my ($cofitfile, $taxonomyId, $taxa, $orth, $maxrank, $mincor) = @_;
    open(COFIT, "<", $cofitfile) || die "Cannot read $cofitfile";
    my $header = <COFIT>;
    chomp $header;
    my @colnames = split /\t/, $header, -1;
    my %colnames = ();
    my $col5 = 0; # File is in 5-column no-header format?
    my $line1 = ""; # if there is no header (to process 1st line in the loop below)
    if ($colnames[0] eq $taxonomyId && scalar(@colnames) == 5) {
        $col5 = 1;
        $line1 = $header;
    } else {
        %colnames = map { $colnames[$_] => $_ } (0..(scalar(@colnames)-1));
        die "Invalid cofitness file $cofitfile"
            unless exists $colnames{"locusId"} && exists $colnames{"hitId"} && exists $colnames{"cofit"} && exists $colnames{"rank"};
    }
    my %cofit = (); # locusId => hitId => [rank,cofitness]
    my %skip = (); # locusIds not in keys of %$orth
    while(my $line = $line1 || <COFIT>) {
	chomp $line;
        $line1 = ""; # do not reprocess
	my @F = split /\t/, $line, -1;
        my ($locusId,$hitId,$cor,$rank);
        if ($col5) {
            die "First field does not match: $F[0]" unless $F[0] eq $taxonomyId;
            (undef,$locusId,$hitId,$rank,$cor) = @F;
        } else {
            die "Error reading cofitness file $cofitfile:\n$line\nhas wrong number of columns"
                unless scalar(@F) == scalar(@colnames);
            $locusId = $F[ $colnames{"locusId"} ];
            $hitId = $F[ $colnames{"hitId"} ];
            $cor = $F[ $colnames{"cofit"} ];
            $rank = $F[ $colnames{"rank"} ];
        }
	die "Colon should not be in locusId in $cofitfile: $locusId" if $locusId =~ m/:/;
	die "Colon should not be in hitId in $cofitfile: $hitId" if $hitId =~ m/:/;

	$locusId = "$taxonomyId:$locusId";
	$hitId = "$taxonomyId:$hitId";
	if (!exists $orth->{$locusId}) {
	    $skip{$locusId} = 1;
	} else {
	    die "Invalid cofitness line $line" unless $rank =~ m/^\d+$/;
	    die "Duplicate hit $hitId for $locusId" if exists $cofit{$locusId}{$hitId};
	    if (exists $orth->{$hitId}) {
                next if defined $maxrank && $rank > $maxrank;
                next if defined $mincor && $cor < $mincor;
		$cofit{$locusId}{$hitId} = [ $rank, $cor ];
	    }
	}
    }
    print STDERR "For taxon $taxonomyId -- read cofitness for " . scalar(keys %cofit) . " with orths and "
        . scalar(keys %skip) . " without\n";
    return \%cofit;
}
