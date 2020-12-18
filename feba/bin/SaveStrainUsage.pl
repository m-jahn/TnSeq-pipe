#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $usage = <<END
Usage: SaveStrainUsage.pl [ -org organism ] [ -fit html/organism ] [ -out g/organism ]

Given the fitness image in the fit directory (from FEBA.R, BarSeqR.pl
or BarSeqtest.pl), saves the strain usage into the g/organism/
directory so that future analyses can use the same set of strains and
genes.
END
    ;

my ($org, $htmldir, $outdir);
die $usage unless GetOptions('org=s' => \$org,
                             'fit=s' => \$htmldir,
                             'out=s' => \$outdir)
    && @ARGV == 0
    && (defined $org || (defined $htmldir && defined $outdir));
$htmldir = "html/$org" unless defined $htmldir;
$outdir = "g/$org" unless defined $outdir;
die "No such directory: $htmldir" unless -d $htmldir;
die "No such directory: $outdir" unless -d $outdir;
die "No such file: $htmldir/fit.image" unless -e "$htmldir/fit.image";

open(R, "|-", "R","--vanilla","--slave") || die "Cannot run R";
print R <<END
  source("$Bin/../lib/FEBA.R");
  load("$htmldir/fit.image");
  SaveStrainUsage(fit, "$outdir");
END
    ;
close(R) || die "Error running SaveStrainUsage() in R";
