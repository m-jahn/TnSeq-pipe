#!/usr/bin/perl -w
# Identify all high fitness values for a set of organisms, regardless of if the experiment succeeded
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable

my $minFit = 2;
my $minT = 5;
my $maxSE = 1;
my $minGMean = 10;
my $maxBelow = 8;

my $usage = <<END
Usage: high_fit.pl -base html org1 ... orgN > high.tab
or
  high_fit.pl -dir dir > high.tab

Optional arguments (with default values shown):
-minFit $minFit -- minimum fitness value to report
-minT $minT -- minimum t value to report
-maxSE $maxSE -- maximum standard error to allow
	(standard error = fit/t)
-maxBelow $maxBelow -- only report fitness values that are
	at least highest_fitness - maxBelow
-minGMean -- only consider experiments where the
	mean reads per gene is at least minGMean.
END
;

my $base = "html";
my $dir;
die $usage
  unless GetOptions('base=s' => \$base,
                    'dir=s' => \$dir,
                    'minFit=f' => \$minFit,
                    'minT=f' => \$minT,
                    'maxSE=f' => \$maxSE,
                    'minGMean=f' => \$minGMean,
                    'maxBelow=f' => \$maxBelow);
die $usage unless (@ARGV > 0) xor (defined $dir);
my @orgs;
my @dirs = ();
if (defined $dir) {
  die "No such directory: $dir\n" unless -d $dir;
  push @dirs, $dir;
} else {
  @orgs = @ARGV;
  foreach my $org (@orgs) {
    die "No such directory: $base/$org\n" unless -d "$base/$org";
    push @dirs, "$base/$org";
  }
}

foreach my $i (0..(scalar(@dirs)-1)) {
  my $dir = $dirs[$i];
  die "No .FEBA.success file in $dir\n" unless -e "$dir/.FEBA.success";
}

print join("\t", qw{orgId locusId sysName desc expName Group Condition_1 Concentration_1 Units_1 Media used gMean maxFit short fit t se})."\n";

foreach my $i (0..(scalar(@dirs)-1)) {
  my $dir = $dirs[$i];
  my $orgId = @orgs > 0 ? $orgs[$i] : "";

  # Load experiment information
  my @exps = ReadTable("$dir/expsUsed", ["name","short","Group","Condition_1","Concentration_1","Units_1","Media"]);
  my @expQ = ReadTable("$dir/fit_quality.tab", ["name","short","gMean","u","maxFit","u"]);
  my %expQ = map { $_->{name} => $_ } @expQ;

  # Load gene information and fitness/t values
  my @genes = ReadTable("$dir/fit_logratios.tab", ["locusId","sysName","desc"]);
  my @tval = ReadTable("$dir/fit_t.tab", ["locusId","sysName","desc"]);

  my $nExp = 0;
  foreach my $exp (@exps) {
    my $expName = $exp->{name};
    next unless exists $expQ{$expName}; # can happen for Time0s with no replicates
    $nExp++;
    my $q = $expQ{$expName};
    my $col = $expName . " " . $exp->{short};
    if (!exists $genes[0]{$col}) {
      print STDERR "Warning: no column for $col in fitness data in $dir\n";
      next;
    }
    die "No such column in t scores: $col" unless exists $tval[0]{$col};

    next unless $q->{gMean} >= $minGMean;
    foreach my $j (0..(scalar(@genes)-1)) {
      my $gene = $genes[$j];
      my $tval = $tval[$j];
      die "log ratios and tvalue tables do not line up"
        unless $gene->{locusId} eq $tval->{locusId};
      my $fit = $gene->{$col};
      next unless $fit >= $minFit && $fit >= $q->{maxFit} - $maxBelow;
      my $t = $tval->{$col};
      next unless $t >= $minT;
      my $se = $fit/$t;
      next unless $se <= $maxSE;
      print join("\t", $orgId, $gene->{locusId}, $gene->{sysName}, $gene->{desc},
                 $expName, $exp->{Group}, $exp->{Condition_1}, $exp->{Concentration_1}, $exp->{Units_1},
                 $exp->{Media}, $q->{u}, $q->{gMean}, $q->{maxFit}, $exp->{short},
                 $fit, $t, $se)."\n";
    }
  }
  print STDERR "Processed $nExp experiments" . ($orgId ? " for $orgId" : "") . "\n";
}

