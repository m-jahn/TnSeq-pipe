#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()
use DBI;

my $metadir = "$Bin/../metadata";
my $gdir = "g";
my $xrefs = "img.xrefs";
my $subsysfile = "subsys.txt";
my $publicationsFile = 'Publications.tab';

my $usage = <<END
Usage: db_setup.pl -db cgi_data/feba.db
        -orth orth_table -orginfo orginfo
        -indir htmldir nickname1 ... nicknameN
Other optional arguments:
    -reanno reannotation_file
    -public public_file
    -metadir $metadir
    -gdir $gdir
    -xrefs $xrefs (from IMGXRefs.pl, or use an empty argument to skip it)
    -subsys $subsysfile
        (SEED subsystems table, or use an empty argument to skip it,
	downloaded from ftp://ftp.theseed.org/subsystems/subsys.txt
	with 4 columns)
    -publications $publicationsFile
    -test -- skip cofitness and strain fitness
    -bad -- include experiments with poor quality metrics in the database

    Sets up the cgi_data/ directory by reading from
    the html directories indir/nickname1 ... indir/nicknameN
    (created by BarSeqR.pl)
    and the g directories gdir/nickname1 ... gdir/nicknameN
    Specifically, creates:
      . the sqlite database (usually cgi_data/feba.db)
      . a BLAST database (usually cgi_data/aaseqs)
      . a usearch database (usually cgi_data/aaseqs.udb)
      . per-strain data files (unless -test is used)
    All of these are put in the same directory as the sqlite database.

    Creates intermediate files db.tablename.org in the working directory.

    The orginfo file should include all of the columns for the Orginfo
    table, but may contain organisms that are not included in the
    command line -- these will not be removed before loading.

    The orth table should include the fields tax1, locus1, tax2, locus2,
    ratio. The loci may be in taxId:locusId format.

    public_file (which is optional) should contain lines of the form
	orgId SetName pubId key1=value1 ... keyN=valueN
	  where the key/value pairs are optional.
	  Use "_" instead of " " to include a space in a value.
      or
        orgId * pubId exps=name1,name2,...,nameN
      For either format, use a pubId of "-" for unpublished
      Other types of (non informative lines):
	comment lines (starting with "#"),
	comments at the end of rules (starting with "#"),
	or empty lines.

    reannotation_file is tab-delimited -- see db_update_reanno.pl for details

    The publications file should be tab-delimited with the fields pubId, title, and URL.
END
;

# global variables for maintaining the work list
my %workRows = (); # table name to number of rows
my $workFile;
my $workTable;
my @workFiles = ();
my @workCommands = ();

sub StartWorkFile($$); # table, file
sub StartWork($$); # table, organism
sub EndWork(); # 
sub WorkPutRow(@); # list of fields to be written as a line
sub WorkPutHash($@); # hash and list of fields to use

# Given a row from the experiment quality table, a row from the experiments table,
# and a list of rules, return the pubId if it matches or undef otherwise
sub ExpToPubId($$$);

{
    my ($dbfile,$indir,$orgfile,$orthfile,$pubfile,$reannofile,$test,$bad);
    (GetOptions('db=s' => \$dbfile,
                'orginfo=s' => \$orgfile,
                'orth=s' => \$orthfile,
                'public=s' => \$pubfile,
                'indir=s' => \$indir,
                'metadir=s' => \$metadir,
                'gdir=s' => \$gdir,
                'xrefs=s' => \$xrefs,
                'reanno=s' => \$reannofile,
                'subsys=s' => \$subsysfile,
                'publications=s' => \$publicationsFile,
                'test' => \$test,
                'bad' => \$bad)
     && defined $indir && defined $orgfile && defined $orthfile && defined $dbfile)
        || die $usage;
    my @orgs = @ARGV;
    die "No such directory: $indir" unless -d $indir;
    die "No such directory: $metadir" unless -d $metadir;
    die "No such directory: $gdir" unless -d $gdir;
    die "No such file: $orgfile" unless -e $orgfile;
    die "No such file: $orthfile" unless -e $orthfile;
    die "No such file: $pubfile" if defined $pubfile && ! -e $pubfile;
    die "No such file: $reannofile" if defined $reannofile && ! -e $reannofile;
    die "No organism nicknames specified\n" if scalar(@orgs) < 1;
    die "No such file: $xrefs" if $xrefs ne "" && ! -e $xrefs;
    if ($subsysfile ne "") {
      die "No such file: $subsysfile" unless -e $subsysfile;
    }
    die "No such file: $publicationsFile" unless -e $publicationsFile;

    # Where to put the db.* files and the BLAST database
    my $outdir = ".";
    if ($dbfile =~ m|/|) {
      $outdir = $dbfile;
      $outdir =~ s!/[^/]+$!!;
    }
    die "Not a directory: $outdir\n" if !-d $outdir;

    # check that all @orgs are unique
    my %orgUniq = map { $_ => 1 } @orgs;
    die "Not all organisms are unique!\n" unless scalar(keys %orgUniq) == scalar(@orgs);

    my $load_kegg = 1;
    my $load_metacyc = 1;
    foreach my $org (@orgs) {
        die "No such directory: $indir/$org" unless -d "$indir/$org";
        foreach my $file (qw{.FEBA.success genes expsUsed fit_quality.tab fit_logratios_good.tab fit_t.tab}) {
            die "Missing file: $indir/$org/$file" unless -e "$indir/$org/$file";
        }
        die "No aaseq2 file for $org in $gdir/$org/aaseq2" unless -e "$gdir/$org/aaseq2";
        die "No genome.fna file for $org in $gdir/$org/aaseq2" unless -e "$gdir/$org/genome.fna";
        die "No swissprot hits for $org in $gdir/blast_results/sprot_$org.m8"
            unless -e "$gdir/blast_results/sprot_$org.m8";
        print STDERR "Warning, no seedanno.tab file in $gdir/$org\n"
            unless -e "$gdir/$org/seedanno.tab";
        if (! -e "$gdir/$org/besthit.kegg") {
            print STDERR "Warning, no kegg besthit file for $org in $gdir/$org/besthit.kegg\n";
            print STDERR "Skipping KEGG\n";
            $load_kegg = 0;
        }
        if (! -e "$gdir/blast_results/metacyc_$org.m8") {
            print STDERR "Warning, no metacyc besthit file for $org in $gdir/$org/metacyc_$org.m8\n";
            print STDERR "Skipping MetaCyc\n";
            $load_metacyc = 0;
        }
    }
    my $formatexe = "$Bin/blast/formatdb";
    die "formatdb not found in $Bin/blast" unless -e $formatexe;

    my @keggtables = qw{ECInfo KEGGCompound KEGGConf KEGGMap};
    foreach my $keggtab (@keggtables) {
        my $file = "$Bin/../kegg/$keggtab";
        die "No such file: $file" unless -e $file;
        push @workCommands, ".import $file $keggtab";
    }
    push @workCommands, "UPDATE KEGGCompound SET formula = NULL WHERE formula = 'NULL';";
    push @workCommands, "UPDATE KEGGCompound SET mass = NULL WHERE mass = 'NULL';";

    print STDERR "Reading " . scalar(@orgs) . " organisms from $indir\n";
    print STDERR "Warning: including bad experiments\n"
      if defined $bad;

    my $tmpdir = $ENV{TMPDIR} || "/tmp";
    my $tmpdbfile = "$tmpdir/db.$$.sqlite3";
    # make sure we can write to the dbfile
    open(DB, ">", $dbfile) || die "Cannot write to $dbfile";
    close(DB) || die "Error writing to $dbfile";
    unlink($dbfile);

    # create the temporary database
    my $sqlfile = "$Bin/../lib/db_setup_tables.sql";
    die "No such file: $sqlfile" unless -e $sqlfile;
    print STDERR "Creating sqlite3 database $tmpdbfile\n";
    unlink($tmpdbfile); # in case it somehow exists
    system("sqlite3 $tmpdbfile < $sqlfile") == 0 || die $!;
    die "Could not create temporary database $tmpdbfile" unless -e $tmpdbfile;

    # Load orginfo
    my @orgfields = qw{name division genus species strain taxonomyId};
    my @orginfo = &ReadTable($orgfile, @orgfields);
    my %orgSeen = map { $_->{name} => 1 } @orginfo;
    my %orgUse = map { $_ => 1 } @orgs;
    foreach my $org (@orgs) {
        die "Organism $org is not described in $orgfile"
            unless exists $orgSeen{$org};
    }

    my @publications_fields = qw{pubId title URL};
    my @publications = &ReadTable($publicationsFile, \@publications_fields);
    my %pubId = map { $_->{pubId} => 1 } @publications;

    my %pubrules = (); # orgId => list of [ rule, pubId ]
    my %ignorePubOrgs = ();
    # each rule is a hash that includes SetName and optionally other fields as well,
    # or, it includes expNames as a list
    if (defined $pubfile) {
        open(PUB, "<", $pubfile) || die "Cannot read $pubfile";
        my $nRules = 0;
        while (my $line = <PUB>) {
            chomp $line;
            $line =~ s/#.*//; # remove comments
            $line =~ s/\s+$//; # remove trailing white space
            $line =~ s/^\s+//; # remove leading white space
            next if $line eq "";
            # rule must have at least orgId and SetName
            my @parts = split /\s+/, $line;
            die "Not enough parts for public rule in\n$line\n" unless scalar(@parts) >= 3;
            my $orgId = shift @parts;
            my $set = shift @parts;
            my $pubId = shift @parts;
            $pubId = "" if $pubId eq "-";
            die "Unknown publication id $pubId" if $pubId ne "" && !exists $pubId{$pubId};
            if (!exists $orgUse{$orgId}) {
              $ignorePubOrgs{$orgId} = 1;
            } else {
              my $rule;
              if ($set eq "*") {
                my $expsSpec = shift @parts;
                die "Invalid public rule $line\n"
                  unless defined $expsSpec && $expsSpec =~ m/^exps=/;
                $expsSpec =~ s/^exps=//;
                my @expNames = split /,/, $expsSpec;
                foreach my $expName (@expNames) {
                  die "Invalid exps= specifier $expsSpec in $pubfile\n"
                    unless $expName =~ m/^[a-zA-Z][a-zA-Z0-9]*\d$/;
                }
                $rule = { 'expNames' => \@expNames };
              } else {
                $rule = { "SetName" => $set };
                foreach my $part (@parts) {
                    my ($key,$value) = split /=/, $part;
                    die "Invalid public assignment $part in\n$line\n"
                        unless defined $value && $value ne "SetName";
                    $rule->{$key} = $value;
                }
              }
              push @{ $pubrules{$orgId} }, [ $rule, $pubId ];
              $nRules++;
            }
        }
        close(PUB);
        print STDERR "Read $nRules relevant rules from $pubfile\n";
        print STDERR "Ignored rules for unknown organisms: " . join(" ", sort keys %ignorePubOrgs)."\n"
          if scalar(keys %ignorePubOrgs) > 0;
    }

    # Create db.Organism
    StartWorkFile("Organism", "db.Organism");
    foreach my $row (@orginfo) {
        $row->{orgId} = $row->{name};
        WorkPutHash($row, qw{orgId division genus species strain taxonomyId})
            if exists $orgUse{$row->{orgId}};
    }
    EndWork();

    # Create db.Gene.*
    foreach my $org (@orgs) {
        my @genes = &ReadTable("$indir/$org/genes",
                               qw{locusId sysName scaffoldId begin end desc name GC});
        die "No genes for $org" unless @genes > 0;
        StartWork("Gene", $org);
        foreach my $row (@genes) {
            $row->{orgId} = $org;
            $row->{gene} = $row->{name};
            $row->{type} = 1 if !exists $row->{type};
            WorkPutHash($row, qw{orgId locusId sysName scaffoldId begin end type strand gene desc GC});
        }
        EndWork();
    }

    # Create db.ScaffoldSeq
    StartWorkFile("ScaffoldSeq", "db.ScaffoldSeq");
    foreach my $orgId (@orgs) {
      my $fnafile = "g/$orgId/genome.fna";
      open(my $fh, "<", $fnafile) || die "Cannot read $fnafile";
      my $state = {};
      while (my ($scaffoldId, $sequence) = ReadFastaEntry($fh, $state)) {
        $scaffoldId =~ s/\s.*//;
        die "Invalid header $scaffoldId in $fnafile" if $scaffoldId eq "";
        WorkPutRow($orgId, $scaffoldId, $sequence);
      }
      close($fh) || die "Error reading $fnafile";
    }
    EndWork();

    # Create db.Ortholog
    # Do not use ReadTable as it may be quite large
    open(ORTH, "<", $orthfile) || die "Cannot read $orthfile";
    my $orthHeader = <ORTH>;
    chomp $orthHeader;
    my @orthHeaderSeen = split /\t/, $orthHeader;
    my @orthHeader = qw{tax1 locus1 tax2 locus2 ratio};
    die "Not enough fields in $orthfile" unless @orthHeaderSeen >= scalar(@orthHeader);
    foreach my $i (0..$#orthHeader) {
        die "Field $orthHeaderSeen[$i] is named $orthHeader[$i] instead"
            unless $orthHeader[$i] eq $orthHeaderSeen[$i];
    }
    my %orgHasOrth = ();
    StartWorkFile("Ortholog", "db.Ortholog");
    while(<ORTH>) {
        chomp;
        my ($tax1,$locus1,$tax2,$locus2,$ratio) = split /\t/, $_, -1;
        die "Cannot parse ortholog file $_" unless defined $ratio && $ratio =~ m/^[0-9.]+$/;
        next unless exists $orgUse{$tax1} && exists $orgUse{$tax2};
        $orgHasOrth{$tax1} = 1;
        # remove orgId: prefix at beginning of locusIds if necessary
        $locus1 =~ s/^.*://;
        $locus2 =~ s/^.*://;
        WorkPutRow($tax1,$locus1,$tax2,$locus2,$ratio);
    }
    close(ORTH) || die "Error reading $orthfile";
    foreach my $org (@orgs) {
        print STDERR "Warning: No orthologs for $org\n"
            unless exists $orgHasOrth{$org};
    }
    EndWork();

    # Create db.Experiment.*
    # Track which experiments are being kept
    my %expKept = (); # orgId => name => 1
    foreach my $org (@orgs) {
        my @q = &ReadTable("$indir/$org/fit_quality.tab",
                           qw{name short t0set num nMapped nPastEnd nGenic nUsed gMed gMedt0 gMean
                                cor12 mad12 mad12c mad12c_t0 opcor adjcor gccor maxFit u});
        # exps has additional fields fields
        my @exps = &ReadTable("$indir/$org/expsUsed",
                              qw{name SetName Date_pool_expt_started Person Mutant.Library Description
                                 Index Media Growth.Method Group
                                 Condition_1 Units_1 Concentration_1
                                 Condition_2 Units_2 Concentration_2});
        my %exps = map {$_->{name} => $_} @exps;
        # fields that are often in exps but are not enforced
        my @optional = qw{MediaStrength Temperature pH Shaking Growth.Method Liquid.v..solid Aerobic_v_Anaerobic
                          Growth.Plate.ID Growth.Plate.wells};
        foreach my $field (@optional) {
            if (!exists $exps[0]{$field}) {
                print STDERR "Field $field is not in $indir/$org/expsUsed -- using blank values\n"
                  unless $field eq "MediaStrength"; # usually missing
                while (my ($name,$exp) = each %exps) {
                    $exp->{$field} = "";
                }
            }
        }

        my $nExpsRead = scalar(@q);
        # only successful experiments
        @q = grep { $_->{u} eq "TRUE" } @q
          unless defined $bad;
        my %q = map { $_->{name} => $_ } @q;
        my $nExpsSucceeded = scalar(@q);
        if (defined $pubfile) {
          foreach my $q (@q) {
            $q->{pubId} = &ExpToPubId($q, $exps{$q->{name}}, $pubrules{$org});
          }
          @q = grep { defined $_->{pubId} } @q;
          # Also warn if experiments are specified to be released but do not exist
          foreach my $rule (@{ $pubrules{$org} }) {
            if (exists $rule->[0]{expNames}) {
              foreach my $expName (@{ $rule->[0]{expNames} }) {
                print STDERR "Warning: org $org experiment $expName is to be released, but does not exist or was low quality!\n"
                  unless exists $q{$expName};
              }
            }
          }
        }
        my $nFiltered = scalar(@q);
        print STDERR "$org: read $nExpsRead experiments, $nExpsSucceeded kept";
        print STDERR ", filtered to $nFiltered" if $nFiltered < $nExpsSucceeded;
        print STDERR "\n";
        print STDERR "Warning, no experiments to show for $org\n" if $nFiltered == 0;

        StartWork("Experiment",$org);
        foreach my $row (@q) {
            $expKept{$org}{$row->{name}} = 1;
            $row->{"orgId"} = $org;
            $row->{"expName"} = $row->{"name"};
            $row->{"expDesc"} = $row->{"short"};
            $row->{"timeZeroSet"} = $row->{"t0set"};

            my $id = $row->{name};
            my $exp = $exps{$id} || die "No matching metadata for experiment $org $id";
            # put fields from exps into output with a new name (new name => original name)
            my %remap = ( "expDescLong" => "Description",
                          "mutantLibrary" => "Mutant.Library",
                          "expGroup" => "Group",
                          "dateStarted" => "Date_pool_expt_started",
                          "setName" => "SetName",
                          "seqindex" => "Index",
                          "temperature" => "Temperature",
                          "shaking" => "Shaking",
                          "vessel" => "Growth.Method",
                          "aerobic" => "Aerobic_v_Anaerobic",
                          "liquid" => "Liquid.v..solid",
                          "pH" => "pH",
                          "growthPlate" => "Growth.Plate.ID",
                          "growthWells" => "Growth.Plate.wells",
                          "nGenerations" => "Total.Generations"
                );
            # and lower case these fields
            foreach my $field (qw{Person Media
                                 Condition_1 Units_1 Concentration_1
                                 Condition_2 Units_2 Concentration_2
                                 Condition_3 Units_3 Concentration_3
                                 Condition_4 Units_4 Concentration_4}) {
                $remap{lc($field)} = $field;
            }
            while (my ($new,$old) = each %remap) {
                $row->{$new} = !defined $exp->{$old} || $exp->{$old} eq "NA" ? "" : $exp->{$old};
            }
            # and compute nGenerations from StartOD and EndOD
            if ($row->{nGenerations} eq ""
                && $exp->{StartOD} ne "" && $exp->{StartOD} ne "NA" && $exp->{StartOD} > 0
                && $exp->{EndOD} ne "" && $exp->{EndOD} ne "NA" && $exp->{EndOD} > 0) {
              my $lr = log( $exp->{EndOD} / $exp->{StartOD} ) / log(2);
              $row->{nGenerations} = $lr if $lr > 0;
            }
            # ad-hoc fixes for issues in the spreadsheets that the metadata was entered in:
            # hidden spaces at end, condition_ = "None" instead of empty, or expGroup not lowercase
            foreach my $i (1..4) {
              my $field = "condition_" . $i;
              $row->{$field} =~ s/ +$//;
              $row->{$field} = "" if $row->{$field} eq "None";
            }
            $row->{expGroup} = lc($row->{expGroup}) unless $row->{expGroup} eq "pH";

            $row->{pubId} = "" if !defined $row->{pubId};

            # Compute mediaStrength from MediaStrength, which should have an x suffix
            my $strength = $exp->{MediaStrength};
            $strength =~ s/^\s+//;
            $strength =~ s/\s+$//;
            $strength = "" if $strength eq "NA";
            if ($strength ne "") {
              if ($strength =~ m/^(\d*[.]?\d*)x$/i && $1 ne "." && $1 ne "") {
                $row->{mediaStrength} = $1;
              } else {
                print STDERR "Warning: cannot parse media strength $exp->{MediaStrength} for experiment $exp->{orgId} $exp->{expName} $exp->{expDesc} -- it should be something like 0.9x or 1x\n";
              }
            }
            $row->{mediaStrength} = 1 if !defined $row->{mediaStrength};
            WorkPutHash($row, qw{orgId expName expDesc timeZeroSet num nMapped nPastEnd nGenic
                                 nUsed gMed gMedt0 gMean cor12 mad12 mad12c mad12c_t0
                                 opcor adjcor gccor maxFit
                                 expGroup expDescLong mutantLibrary person dateStarted setName seqindex
                                 media mediaStrength
                                 temperature pH vessel aerobic liquid shaking
                                 condition_1 units_1 concentration_1
                                 condition_2 units_2 concentration_2
                                 condition_3 units_3 concentration_3
                                 condition_4 units_4 concentration_4
                                 growthPlate growthWells nGenerations pubId});
        }
        EndWork();
    }

    # Create db.Publications
    StartWorkFile("Publication", "db.Publication");
    foreach my $pub (@publications) {
      WorkPutHash($pub, qw{pubId title URL});
    }
    EndWork();

    # Create db.GeneFitness.* and, if enough experiments, db.Cofit.*
    my %orgCofit = (); # org => 1 if organism has cofitness
    foreach my $org (@orgs) {
        my $fit_file = defined $bad ? "$indir/$org/fit_logratios.tab"
          : "$indir/$org/fit_logratios_good.tab";
        my @fitNames = &ReadColumnNames($fit_file);
        my @colNamesFull = grep m/^set/i, @fitNames; # must start with set or Set
        print STDERR "Warning, no successful experiments for $org\n"
          if @colNamesFull == 0;
        my @colNames = @colNamesFull;
        map { s/ .*$//; } @colNames; # skip annotations after the name
        my @fit = &ReadTable($fit_file, "locusId");
        my @t = &ReadTable("$indir/$org/fit_t.tab",grep {$_ ne "comb"} @fitNames);

        StartWork("GeneFitness",$org);
        my %geneFit = (); # geneId => list
        my $expKept = $expKept{$org};
        foreach my $iRow (0..$#fit) {
            my $fit_row = $fit[$iRow];
            my $t_row = $t[$iRow];
            my $locusId = $fit_row->{locusId};
            die "Mismatched locusIds in $fit_file and fit_t.tab, row $iRow"
                unless $locusId eq $t_row->{locusId};
            my @fit_vec = ();
            foreach my $iCol (0..$#colNamesFull) {
                my $colFull = $colNamesFull[$iCol];
                die "Missing from fit: $iCol $colFull" unless defined $fit_row->{$colFull};
                die "Missing from t: $iCol $colFull" unless defined $t_row->{$colFull};
                die "Missing from colNames: $iCol $colFull" unless defined $colNames[$iCol];
                # censor this column from the database and the fitness vectors used for cofitness
                next unless exists $expKept->{ $colNames[$iCol] };
                WorkPutRow($org, $locusId, $colNames[$iCol],
                           $fit_row->{$colFull}, $t_row->{$colFull});
                push @fit_vec, $fit_row->{$colFull};
            }
            $geneFit{$locusId} = \@fit_vec;
        }
        EndWork();

        my @tophits = (); # a list of rows
        my @colsKept = grep exists $expKept->{$_}, @colNames;
        my $nExpKept = scalar(@colsKept);
        if (defined $test) {
          print STDERR "Skipping cofitness -- testing\n";
        } elsif ($nExpKept < 15) {
            print STDERR "Not computing cofitness for $org -- just $nExpKept experiments\n";
        } else {
            print STDERR "Computing top hits for $org\n";
            my $fit_file = "db.fittab.$org";
            my $cofit_file = "db.cofit_tab.$org";
            $orgCofit{$org} = 1;

            open(FIT, ">", $fit_file) || die "Error writing to $fit_file";
            print FIT join("\t","locusId",@colsKept)."\n";
            # go through in a reasonable order in case it helps indexing the db
            # (results from TopCofit.R will be in this same order)
            foreach my $locusId (map { $_->{locusId} } @fit) {
                print FIT join("\t", $locusId, @{ $geneFit{$locusId} })."\n";
            }
            close(FIT) || die "Error writing to $fit_file";

            system("$Bin/TopCofit.R", $fit_file, $cofit_file) == 0
                || die "Error running $Bin/TopCofit.R $fit_file $cofit_file\nStatus: $!";
            @tophits = &ReadTable($cofit_file, qw(locusId hitId rank cofit));
            unlink($fit_file);
            unlink($cofit_file);
        }
        StartWork("Cofit", $org);
        foreach my $row (@tophits) {
            $row->{orgId} = $org;
            WorkPutHash($row, qw{orgId locusId hitId rank cofit});
        }
        EndWork();
    } # end do GeneFitness and Cofit for each organism

    # Create db.SpecificPhenotype.*
    foreach my $org (@orgs) {
        my $specfile = "$indir/$org/specific_phenotypes";
        if (! -e $specfile) {
            print STDERR "No specific_phenotypes files for $org (usually means too few experiments)\n";
            next;
        }
        my @spec = &ReadTable($specfile, qw{locusId name});
        StartWork("SpecificPhenotype",$org);
        foreach my $row (@spec) {
            $row->{orgId} = $org;
            $row->{expName} = $row->{name};
            WorkPutHash($row, qw{orgId expName locusId})
                if exists $expKept{$org}{ $row->{expName} };
        }
        EndWork();
    }

    # Create db.GeneDomain.*
    foreach my $org (@orgs) {
        # my @pFam = &ReadTable("$indir/g/$org/pfam.tab",
        my @pFam = &ReadTable("g/$org/pfam.tab",
                              qw{locusId domainId domainName begin end score evalue});
        # my @tigrFam = &ReadTable("$indir/g/$org/tigrfam.tab",
        my @tigrFam = &ReadTable("g/$org/tigrfam.tab",
                                 qw{locusId domainId domainName begin end score evalue});

        die "No Fam data for $org" unless @pFam > 0 or @tigrFam > 0;

        
        # fields that are usually in tigrfam but are not enforced
        # my @tigrInfo = &ReadTable("$indir/tigrinfo", qw{id tigrId type roleId geneSymbol ec definition});
        my @tigrInfo = &ReadTable("tigrinfo", qw{id tigrId type geneSymbol ec definition});
        my %tigrInfo = map {$_->{tigrId} => $_} @tigrInfo;

        StartWork("GeneDomain",$org);
        foreach my $row (@pFam) {
            $row->{"domainDb"} = 'PFam';
            $row->{"orgId"} = $org;
            $row->{"domainId"} =~ s/\.\d+//; # remove the suffix to allow easier searches
            $row->{"type"} = "";
            $row->{"geneSymbol"} = "";
            $row->{"ec"} = "";
            $row->{"definition"} = "";
            WorkPutHash($row, qw{domainDb orgId locusId domainId domainName begin end score evalue type geneSymbol ec definition});
        };

        foreach my $row (@tigrFam) {
            $row->{"domainDb"} = 'TIGRFam';
            $row->{"orgId"} = $org;

            my $id = $row->{domainId};
            my $info = $tigrInfo{$id};
            print "No matching metadata for tigrFam $org $id $row->{locusId}\n" if !defined $info;

            # put fields from info into the output
            if (defined $info) {
                # $row->{"locusId"} = $info->{"id"};
                # $row->{"domainId"} = $info->{"tigrId"};
                $row->{"type"} = $info->{"type"};
                $row->{"geneSymbol"} = $info->{"geneSymbol"};
                $row->{"ec"} = $info->{"ec"};
                $row->{"definition"} = $info->{"definition"};
            } else {
                $row->{"type"} = "";
                $row->{"geneSymbol"} = "";
                $row->{"ec"} = "";
                $row->{"definition"} = "";
            }
            WorkPutHash($row, qw{domainDb orgId locusId domainId domainName begin end score evalue type geneSymbol ec definition});
        }
        EndWork();
    }

    # Strain fitness tables are expected to be sorted by scaffoldId and position (but this is not checked).
    # Strains with scaffold="pastEnd" or enoughT0 != "TRUE" are ignored.
    # Strain fitness values are rounded to just 1 decimal point.
    # Strain fitness tables are indexed by seek positions to the kb, which is in StrainDataSeek table.
    if (defined $test) {
      print STDERR "Skipping strain fitness, testing\n";
    } else {
      foreach my $org (@orgs) {
        open(OUT, ">", "db.StrainFitness.$org") || die "Cannot write to db.StrainFitness.$org";
        my $infile = "$indir/$org/strain_fit.tab";
        if (!open(IN, "<", $infile)) {
          print STDERR "Cannot read $infile -- skipping strain fitness for $org\n";
        } else {
          StartWork("StrainDataSeek", $org);
          my $headerLine = <IN>;
          chomp $headerLine;
          my @colNames = split /\t/, $headerLine;
          # column name => index
          my %colNames = map {$colNames[$_] => $_} (0..(scalar(@colNames)-1));
          my @metaShow = qw{barcode scaffold strand pos locusId used};
          foreach my $col (@metaShow) {
            die "No column for $col in $infile" unless exists $colNames{$col};
          }
          foreach my $col (qw{enoughT0}) {
            die "No column for $col in $infile" unless exists $colNames{$col};
          }
          my $expKept = $expKept{$org};
          foreach my $expName (keys %$expKept) {
            die "No column for $expName in $infile" unless exists $colNames{$expName};
          }
          my @exps = grep { exists $expKept->{$_} } @colNames;
          my @expI = map { $colNames{$_} } @exps;
          
          print OUT join("\t", @metaShow, @exps)."\n";
          my @metaShowI = map { $colNames{$_} } @metaShow;
          my $iColEnoughT0 = $colNames{enoughT0};
          my $iColScaffold = $colNames{scaffold};
          my $iColLocus = $colNames{locusId};
          my $iColPos = $colNames{pos};
          
          my $lastScaffold = "";
          my $lastKb = -1;
          
          # No hash look ups in inner loop make this ~2x faster.
          while (<IN>) {
            chomp;
            my @F = split /\t/, $_, -1;
            die "Wrong number of rows in\n$_\nfrom $infile" unless scalar(@F) == scalar(@colNames);
            my $scaffold = $F[$iColScaffold];
            next unless $F[$iColEnoughT0] eq "TRUE" && $scaffold ne "pastEnd";
            $F[$iColLocus] = "" if $F[$iColLocus] eq "NA";
            my @out = @F[ @metaShowI ];
            push @out, map { sprintf("%.1f",$_) } @F[@expI];
            
            my $kb = int($F[$iColPos] / 1000.0);
            if ($scaffold ne $lastScaffold || $kb != $lastKb) {
              WorkPutRow($org, $scaffold, $kb, tell(OUT));
            }
            print OUT join("\t", @out)."\n";
            
            $lastScaffold = $scaffold;
            $lastKb = $kb;
          }
          EndWork();
          close(IN) || die "Error reading $infile";
        } # end if have strain data for org
        close(OUT) || die "Error writing to db.StrainFitness.$org";
        print STDERR "Wrote db.StrainFitness.$org\n";
      } # end loop over orgs
    } # end else (not test mode)

    # Load SwissProt hits
    if (defined $outdir) {
        my $bhfile = "$outdir/db.BestHitSwissProt";
        my $swfile = "$outdir/db.SwissProtDesc";
        # do they exist for each organism?
        foreach my $org (@orgs) {
          my $hitsFile = "$gdir/blast_results/sprot_$org.m8";
          die "No such file: $hitsFile" unless -e $hitsFile;
        }
        print STDERR "Running SprotBestHit.pl\n";
        system("$Bin/SprotBestHit.pl","-gdir",$gdir,"-out",$outdir,@orgs) == 0
          || die "SprotBestHit.pl failed";
        die if ! -e $bhfile && -e $swfile;
        push @workCommands, ".import $bhfile BestHitSwissProt";
        push @workCommands, ".import $swfile SwissProtDesc";
    }

    # rename the strain fitness files
    if (defined $outdir && $outdir ne "." && !defined $test) {
      print STDERR "Removing any preexisting db.StrainFitness.* files from $outdir\n";
      system("rm -f $outdir/db.StrainFitness.* >& /dev/null");
      foreach my $org (@orgs) {
        system("mv", "db.StrainFitness.$org", $outdir) == 0 || die "mv db.StrainFitness.$org $outdir failed: $!";
        die "No such file: $outdir/db.StrainFitness.$org" unless -e "$outdir/db.StrainFitness.$org";
      }
      print STDERR "Moved db.StrainFitness.* into $outdir\n";
    }

    push @workCommands, ".import $xrefs LocusXref" if $xrefs ne "";

    # Build the ConservedCofit table
    if (scalar(keys %orgCofit) >= 2) {
        my @command = ("$Bin/conserved_cofit.pl",
                       "-orth", $orthfile,
                       "-out", "db.ConservedCofit",
                       "-table",
                       "-rank", 10,
                       "-cor", 0.6);
        foreach my $org (sort keys %orgCofit) {
            push @command, $org;
            push @command, "db.Cofit.$org";
        }
        system(@command) == 0
            || die "conserved_cofit.pl failed: " . join(" ", @command) . "\n";
        push @workCommands, ".import db.ConservedCofit.pairs ConservedCofit"; 
    }

    # Build the SEED tables
    system("$Bin/db_setup_SEED.pl","-gdir",$gdir,"-subsys",$subsysfile,@orgs) == 0
        || die "$Bin/db_setup_SEED.pl failed";
    foreach my $table (qw{SEEDAnnotation SEEDClass SEEDRoles SEEDAnnotationToRoles}) {
      push @workCommands, ".import db.$table $table";
    }

    # load the other data into sqlite3
    # Run the commands
    open(SQLITE, "|-", "sqlite3", "$tmpdbfile") || die "Cannot run sqlite3 on $tmpdbfile";
    print SQLITE ".bail on\n";
    print SQLITE ".mode tabs\n";
    print STDERR "Loading tables\n";
    foreach my $workCommand (@workCommands) {
      print SQLITE "$workCommand\n";
    }
    close(SQLITE) || die "Error running sqlite3 commands\n";

    # Check #rows in each table
    my $dbh = DBI->connect("dbi:SQLite:dbname=$tmpdbfile", "", "", { RaiseError => 1 }) || die $DBI::errstr;
    while (my ($table, $nRowsExpect) = each %workRows) {
      my ($nRowsActual) = $dbh->selectrow_array("SELECT COUNT(*) FROM $table");
      die "counting rows in $table failed" unless defined $nRowsActual;
      die "Failed to load $table: expect $nRowsExpect rows but see $nRowsActual rows instead\n"
        unless $nRowsActual == $nRowsExpect;
    }

    # Sanity check orthologs:
    my ($nOrthRows) = $dbh->selectrow_array("SELECT COUNT(*) FROM Ortholog");
    my ($nOrthRows1) = $dbh->selectrow_array("SELECT COUNT(*) FROM Ortholog JOIN Gene ON orgId1=orgId AND locusId1=locusId");
    my ($nOrthRows2) = $dbh->selectrow_array("SELECT COUNT(*) FROM Ortholog JOIN Gene ON orgId2=orgId AND locusId2=locusId");
    die "Cannot count orths" unless defined $nOrthRows && defined $nOrthRows1 && defined $nOrthRows2;
    die "Invalid values of locusId1 in Ortholog table: $nOrthRows != $nOrthRows1" unless $nOrthRows == $nOrthRows1;
    die "Invalid values of locusId2 in Ortholog table: $nOrthRows != $nOrthRows1" unless $nOrthRows == $nOrthRows2;

    # Create the FitByExp_org table for each organism
    foreach my $org (@orgs) {
      my $tab = "'" . "FitByExp_" . $org . "'";
      my $create_statement = $dbh->prepare(qq{CREATE TABLE $tab
                                                (expName TEXT NOT NULL,
                                                 locusId TEXT NOT NULL,
                                                 fit REAL NOT NULL,
                                                 t REAL NOT NULL,
                                                 PRIMARY KEY (expName,locusId)); });
      $create_statement->execute() || die "Cannot create $tab";
      my $insert_statement = $dbh->prepare(qq{INSERT INTO $tab SELECT expName,locusId,fit,t
                                                FROM GeneFitness where orgId = ? ORDER BY expName; });
      $insert_statement->execute($org) || die "Cannot fill data into $tab";
      print STDERR "Filled $tab\n";
    }

    $dbh->disconnect();

    # Build the SpecOG table
    system("$Bin/db_setup_specOG.pl", "-db", $tmpdbfile, "-dir", $outdir) == 0
      || die "db_setup_specOG.pl failed";

    system("mv",$tmpdbfile,$dbfile) == 0 || die "mv $tmpdbfile $dbfile failed: $!";
    die "mv $tmpdbfile $dbfile failed" unless -e $dbfile;
    print STDERR "Successfully created $dbfile\n";
    print STDERR "Cleaning up\n";
    # Delete the files
    foreach my $file (@workFiles) {
      unlink($file);
    }

    if ($load_kegg) {
      my @inputs = map "$gdir/$_/besthit.kegg", @orgs;
      system("$Bin/db_setup_kegg.pl", "-db", $dbfile, "-dir", $outdir, @inputs) == 0
        || die "db_setup_kegg.pl failed";
    } else {
      print STDERR "Skipping kegg hits\n";
    }

    if ($load_metacyc) {
      system("$Bin/db_setup_metacyc.pl", "-db", $dbfile, "-dir", $outdir, @orgs) == 0
        || die "db_setup_metacyc.pl failed";
    }

    # Load the media information
    my @mediacmd = ("$Bin/make_media_table.pl", "-metadir", $metadir, "-out", $outdir);
    if (defined $dbfile) {
      push @mediacmd, ("-db", $dbfile);
    } elsif (defined $outdir) {
      push @mediacmd, ("-out", $outdir);
    }
    system(@mediacmd) == 0 || die "Error running\n" . join(" ",@mediacmd) . "\n: $!";

    # Load the reannotation information
    if (defined $reannofile) {
      if (defined $dbfile) {
        my @reannocmd = ("$Bin/db_update_reanno.pl", "-nometacyc", "-reanno", $reannofile, "-db", $dbfile);
        system(@reannocmd) == 0 || die "Error running\n" . join(" ", @reannocmd) . "\n: $!";
      } else {
        print STDERR "Ignoring -reanno $reannofile because no database specified.\n";
      }
    }

    # Compute MetaCyc pathway coverage
    my @covercmd = ("$Bin/db_update_metacyc_coverage.pl", "-db", $dbfile);
    system(@covercmd) == 0 || die "Error running " . join(" ", @covercmd) . "\n: $!";

    # Make the BLAST database
    my @files = map "$gdir/$_/aaseq2", @orgs;
    my $blastdb = "$outdir/aaseqs";
    my $catcmd = "cat " . join(" ",@files) . " > $blastdb";
    system($catcmd) == 0 || die "Error running\n$catcmd\n: $!";
    print STDERR "Formatting $blastdb\n";
    my @formatcmd = ($formatexe, "-p", "T", "-o", "T", "-i", $blastdb);
    system(@formatcmd) == 0 || die "Error running\n".join(" ",@formatcmd)."\n: $!";

    # and make the udb database
    my $udbfile = "$outdir/aaseqs.udb";
    my $usearch = "$Bin/usearch";
    if (-x $usearch) {
      system("$usearch -makeudb_ublast $blastdb -output $udbfile") == 0
        || die "Failed to make the udb file";
    } else {
      print STDERR "Warning: skipped making aaseqs.udb:\n$usearch is not an executable\n";
    }
    print STDERR "Success\n";
}

sub StartWorkFile($$) {
    my ($table,$file) = @_;
    die "Already working on $workTable $workFile when calling StartWorkFile with $table $file"
        if defined $workTable || defined $workFile;
    $workTable = $table;
    $workFile = $file;
    push @workFiles, $workFile;
    push @workCommands, ".import $workFile $workTable";
    open(WORK, ">", $workFile) || die "Cannot write to $workFile";
    $workRows{$workTable} = 0 if !exists $workRows{$workTable};
}

sub StartWork($$) {
    my ($table,$org) = @_;
    StartWorkFile($table,"db.$table.$org");
}

sub EndWork() {
    die "Illegal call to EndWork()" unless defined $workFile && defined $workTable;
    close(WORK) || die "Error writing to $workFile";
    print STDERR "Wrote $workFile\n";
    $workFile = undef;
    $workTable = undef;
}

sub WorkPutRow(@) {
    die "Illegal call to WorkPutRow" unless defined $workFile && defined $workTable;
    print WORK join("\t", @_)."\n" || die "Error writing to $workFile";
    $workRows{$workTable}++;
}

sub WorkPutHash($@) {
    my ($row,@fields) = @_;
    foreach my $field (@fields) {
        die "No such field $field" unless exists $row->{$field};
        die "Undefined field $field" if !defined $row->{$field};
    }
    WorkPutRow(map $row->{$_}, @fields);
}

# Convert white space to underscores
sub wsu($) {
  my ($value) = @_;
  $value =~ s/ /_/g;
  return $value;
}

sub ExpToPubId($$$) {
    my ($q, $exp, $rules) = @_;
    foreach my $row (@$rules) {
      my ($rule, $pubId) = @$row;
      if (exists $rule->{expNames}) {
        foreach my $expName (@{ $rule->{expNames} }) {
          return $pubId if $expName eq $exp->{name};
        }
      } else {
        next unless $rule->{SetName} eq $exp->{SetName};
        # See if all other fields match
        my $match = 1;
        while (my ($key,$value) = each %$rule) {
          if (exists $q->{$key}) {
            $match = 0 unless &wsu($q->{$key}) eq &wsu($value);
          } elsif (exists $exp->{$key}) {
            $match = 0 unless &wsu($exp->{$key}) eq &wsu($value);
          } else {
            die "Cannot handle rule for $key=$value because $key does not exist for either q or experiment for $exp->{name}";
          }
        }
        return $pubId if $match;
      }
    }
    return undef;
}
