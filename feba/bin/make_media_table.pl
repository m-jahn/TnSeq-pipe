#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils;
use Compounds;

my $metadir = "$Bin/../metadata";
my $formatDoc = Compounds::GetCompoundFormatDoc();

my $usage = <<END
make_media_table.pl [ -out . ] [ -db feba.db ]
             [ -metadir $metadir ]
        The metadata driectory should include a media file, a mixes
	file, and a Compounds.tsv file.
        writes db.Compounds and db.MediaComponents to output directory.
        If -db is specified, also loads these into the sqlite3 database.

$formatDoc
END
    ;

{
    my $db;
    my $outdir = ".";

    die $usage unless GetOptions('metadir=s' => \$metadir,
                                 'out=s' => \$outdir,
                                 'db=s' => \$db)
        && @ARGV == 0;
    die "Not a directory: $metadir" unless -d $metadir;
    die "Not a directory: $outdir" unless -d $outdir;
    die "No such file: $db" if defined $db && ! -e $db;

    LoadCompounds($metadir);
    LoadMedia($metadir);

    my @unknownComponents = GetUnknownComponents();
    my $nNoSyn = scalar(@unknownComponents);
    if ($nNoSyn > 0) {
        print STDERR "Unrecognized compounds: $nNoSyn\n";
        print STDERR join("\t", @unknownComponents)."\n";
    }

    my @undefMedia = GetUndefMedia();
    my $nUndefMedia = scalar(@undefMedia);
    if ($nUndefMedia > 0) {
        print STDERR "Media with no definitions: $nUndefMedia\n";
        print STDERR join("\t", sort @undefMedia)."\n";
    }

    WarnReusedComponents();

    # write out the compounds table
    my $DbCompoundsFile = "$outdir/db.Compounds";
    open(COMPOUNDS, ">", $DbCompoundsFile) || die "Cannot write to $DbCompoundsFile";
    foreach my $compound (GetCompoundList()) {
        print COMPOUNDS join("\t", $compound, GetCompoundMW($compound), GetCompoundCAS($compound))."\n";
    }
    close(COMPOUNDS) || die "Error writing to $DbCompoundsFile";
    print STDERR "Wrote $DbCompoundsFile\n";

    # write out all the media definitions
    my $DbComponentsFile = "$outdir/db.MediaComponents";
    open(COMPONENTS, ">", $DbComponentsFile) || die "Cannot write to $DbComponentsFile";
    foreach my $media (GetMedias) {
        my $list = GetMediaComponents($media);
        foreach my $row (@$list) {
            my ($compound,$concentration,$units,$mix) = @$row;
            # database schema allows for concentration or units to be missing, but this script does not
            print COMPONENTS join("\t", $media, $compound, $concentration, $units, $mix)."\n";
        }
    }
    foreach my $mix (GetMixes) {
      my $list = GetMixComponents($mix);
      foreach my $row (@$list) {
        my ($compound, $concentration, $units) = @$row;
        my $mix2 = "";
        print COMPONENTS join("\t", $mix, $compound, $concentration, $units, $mix2)."\n";
      }
    }

    close(COMPONENTS) || die "Error writing to $DbComponentsFile";
    print STDERR "Wrote $DbComponentsFile\n";

    if (defined $db) {
        print STDERR "Loading tables into $db\n";
        my @commands = (".bail on",
                        ".mode tabs",
                        "DELETE FROM Compounds;",
                        ".import $DbCompoundsFile Compounds",
                        "DELETE FROM MediaComponents;",
                        ".import $DbComponentsFile MediaComponents",
                        ".headers on",
                        "SELECT count(distinct media) nMedia, count(*) nTotalComponents from MediaComponents;",
                        "SELECT count(*) nCompounds from Compounds;"
            );
        open(SQLITE, "|-", "sqlite3", "$db") || die "Cannot run sqlite3 on $db";
        foreach my $command (@commands) {
            print SQLITE "$command\n";
        }
        close(SQLITE) || die "Error running sqlite3";
    }
}
