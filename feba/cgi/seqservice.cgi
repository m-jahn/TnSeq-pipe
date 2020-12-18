#!/usr/bin/perl -w
#######################################################
## seqservice.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: seq -- the protein sequence to look for matches to.
# Optional parameters: maxhits -- default is 20. Cannot be raised above 50.
# debug -- write lines starting with # to output to record status
# html -- return HTML code (with no wrappers) for the result, instead of a table
# wrap -- return HTML code with HTML tags (for testing). Must specify with html.
#
# Must start the ublast service first with bin/start_ublast_service.pl
#
# By default, returns a tab-delimited table of FastBLAST hits, or, a
# simple two-line table with Error and the error (or "no hits").
#
# For an example of a page that uses this service see
# ../images/fitblast_example.html
# ../images/fitblast.js

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use Cwd; # for absolute paths
use Time::HiRes qw(usleep);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
if ($cgi->param('wrap') && $cgi->param('html')) {
    print $cgi->header('text/html');
} else {
    print $cgi->header(-type => 'text/plain',  
                       # allow cross-site scripts to use this
                       -access_control_allow_origin => '*');
}

my $seq = $cgi->param('seq');
$seq =~ s/[ \r\n]//g;
if (! $seq) {
    print "Error\nNo sequence specified\n";
    exit(0);
}
# Allow * as internal stop codon
unless ($seq =~ m/^[A-Z*]+$/) {
    print "Error\nIllegal sequence\n";
    exit(0);
}

my $maxhits = $cgi->param('maxhits');
if (defined $maxhits && $maxhits > 0) {
    $maxhits = $maxhits + 0; # convert to integer
} else {
    $maxhits = 20;
}

my $debug = $cgi->param('debug');

my $base = "../cgi_data/ublast_service";
die "No such directory: $base" unless -d $base;

my $len = length($seq);
my $dir = "$base/seq$$.$len";
mkdir($dir);
die "Cannot mkdir $dir" unless -d $dir;
system("chmod a+w $dir");

open(FASTA, ">", "$dir/faa") || die "Cannot write to $dir/faa";
print FASTA ">seq\n$seq\n";
close(FASTA) || die "Error writing to $dir/faa";

my $path = getcwd();
# note $dir.q is a file in the main ublast_service directory
open(QFILE, ">", "$dir.q") || die "Cannot write to $dir.q";
print QFILE "$path/$dir/faa $path/$dir\n";
close(QFILE) || die "Error writing to $dir.q";

# poll every 0.05 seconds for up to 15 seconds
for (my $i = 0; $i < 300; $i++) {
    if (-e "$dir/done.q") {
	last;
    }  else {
	usleep(50*1000); # in microseconds
    }
}
if (! -e "$dir/done.q") {
    print "Error\nTimeout\n";
    exit;
}

#else
my $dbh = Utils::get_dbh() || die "Cannot access database";
my @rows = ();
open(OUT, "<", "$dir/ublast.out") || die "No $dir/ublast.out file";
while(<OUT>) {
    chomp;
    my @F = split /\t/, $_;
    push @rows, \@F;
    last if @rows >= $maxhits;
}
close(OUT);

unlink("$dir/ublast.out");
unlink("$dir.x");
unlink("$dir.q");
unlink("$dir/faa");
unlink("$dir/done.q");
rmdir("$dir");

my $base_url = url();
$base_url =~ s!/[a-zA-Z_.]+$!!;
my $blastURL = "$base_url/mySeqSearch.cgi?query=$seq";

my $html = $cgi->param('html');
# d3.tsv treats an empty table as an error, so instead, return an actual useful error code
if (@rows == 0 ) {
    if ($html) {
        print qq{No close hits (<A HREF="$blastURL" TITLE="View more hits">more</A>)\n};
    } else {
        print "Error\nNo hits\n";
    }
    exit(0);
}


print join("\t", qw{orgId organism locusId sysName name description identity coverage evalue bits minfit maxfit minT maxT maxcofit})."\n" unless $html;
my $orginfo = Utils::orginfo($dbh);
my @save = ();
foreach my $row (@rows) {
    my ($query,$locusspec,$identity,$alnlen,$mm,$gap,$qBeg,$qEnd,$lBeg,$lEnd,$evalue,$bits) = @$row;
    my ($orgId,$locusId) = split /:/, $locusspec;
    die "Invalid locus $locusspec" unless defined $orgId && defined $locusId;
    my ($sysName,$name,$desc) = $dbh->selectrow_array("SELECT sysName, gene, desc FROM Gene WHERE orgId = ? AND locusId = ?",
					  {}, $orgId, $locusId);
    $sysName = "" if !defined $sysName;
    $name = "" if !defined $name;
    my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(qq{ SELECT min(fit), max(fit), min(t), max(t)
                                                                 FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
							     {}, $orgId, $locusId);
    $minFit = "" unless defined $minFit;
    $maxFit = "" unless defined $maxFit;
    $minT = "" unless defined $minT;
    $maxT = "" unless defined $maxT;
    my $maxCofit = "";
    if (defined $minFit) {
	$maxCofit = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
					  {}, $orgId, $locusId);
        $maxCofit = "" if !defined $maxCofit; # possible if there is no cofitness for this bug
    }
    my $coverage = $alnlen/length($seq);
    if ($html) {
        push @save, [ $orgId, $orginfo->{$orgId}{genome}, $locusId, $sysName, $name, $desc,
                      $identity, $coverage,
                      $minFit, $maxFit, $minT, $maxT, $maxCofit ];
    } else {
        print join("\t",$orgId, $orginfo->{$orgId}{genome},
                   $locusId, $sysName, $name, $desc,
                   $identity,$coverage,$evalue,$bits,
                   $minFit,$maxFit,$minT,$maxT,$maxCofit)."\n";
    }
}

$dbh->disconnect();

if ($html) {
    # mimic fitblast_short() from fitblast.js
    my $mincoverage = 0.75; # hits below this coverage are ignored
    my $mincloseid = 80; # %identity to be considered a close hit
    my $minabsstrong = 2.0;
    # a hit is considered likely to be useful if either
    # cofitness is above threshold or (both absfit and absT are above thresholds)
    my $mincofit = 0.75;
    my $minabsfit = 1.0;
    my $minabsT = 4.0;
    my @pieces = (); # pieces of HTML, 1 for item shown so far

    foreach my $row (@save) {
        my ($orgId, $orgName, $locusId, $sysName, $name, $desc,
            $identity, $coverage, $minFit, $maxFit, $minT, $maxT, $maxCofit) = @$row;
        next unless $coverage >= $mincoverage;
        my $close = $identity >= $mincloseid;
        my $hascofit = $maxCofit >= $mincofit;
        my $hassig = ($minFit <= -$minabsfit && $minT <= $minabsT)
            || ($maxFit >= $minabsfit && $maxT >= $minabsT);
        my $hasstrong = $hassig && ($minFit <= -$minabsstrong || $maxFit >= $minabsstrong);
        my $useful = $hascofit || $hassig;
        if (($close && scalar(@pieces) == 0) || $useful) {
            my $sigstring = $hasstrong ? "strong phenotype" :
                ($hassig ? "has phenotype" : ($minFit ne "" ? "has data" : "no data"));
            if ($minFit ne "") {
                $sigstring = sprintf(qq{<A TITLE="fitness %.1f to %.1f" HREF="$base_url/singleFit.cgi?orgId=%s&locusId=%s">%s</A>},
                                     $minFit, $maxFit, $orgId, $locusId, $sigstring);
            }
            if ($hascofit) {
                $sigstring = sprintf(qq{%s, <A TITLE="top cofitness %.2f" HREF="$base_url/cofit.cgi?orgId=%s&locusId=%s">cofit</A>},
                                     $sigstring, $maxCofit, $orgId, $locusId);
            }
            push @pieces, sprintf(qq{%d%% id. to <A HREF="$base_url/singleFit.cgi?orgId=%s&locusId=%s">%s</A> from %s: %s\n},
                                  int($identity+0.5), $orgId, $locusId,
                                  $name eq "" ? ($sysName eq "" ? $locusId : $sysName) : $name,
                                  $orgName, $sigstring);
        }
        last if $useful;
    }
    print join("<BR>",@pieces) 
        . qq{ (<A HREF="$blastURL" TITLE="view all hits">more</A>)} . "\n";
}


exit(0);

