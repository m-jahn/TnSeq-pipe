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
# Starts the usearch server that supports seqservice.cgi

use strict;

my $debug = 0;
if (@ARGV > 0 && $ARGV[0] eq "-debug") {
    $debug = 1;
    shift @ARGV;
    print STDERR "Debug mode\n";
}
die "Usage: start_ublast_service.pl [ -debug ] cgi_data_directory &\n"  unless @ARGV==1 && -d $ARGV[0];
my $datadir = $ARGV[0];
die "No aaseqs file in $datadir" unless -e "$datadir/aaseqs";
my $usearch = "$datadir/../bin/usearch6";
die "No such exectable: $usearch" unless -x $usearch;
my $basedir = "$datadir/ublast_service";
system("rm","-Rf",$basedir);
mkdir($basedir);
die "Could not make $datadir/ublast_service" unless -d $basedir;
system("chmod a+w $basedir");
system("sed -e 's/[*]/X/g' < $datadir/aaseqs > $basedir/aaseqs");

my $pids = `pgrep usearch6`;
chomp $pids;
print STDERR "Warning, usearch6 is already running: $pids\n" if $pids;

unless($debug) {
    print STDERR "Redirecting to silence usearch6 (use -debug to turn this off)\n";
    print STDERR "usearch will watch $basedir for *.q files\n";
    open(STDIN, "<", "/dev/null");
    open(STDOUT, ">", "/dev/null");
    open(STDERR, ">", "/dev/null");
}

exec($usearch, "-ublast", $basedir, "-server",
     "-db", "$basedir/aaseqs",
     "-query_seqtype","amino",
     "-maxhits",50,"-maxaccepts",50,
     "-blast6out","ublast.out","-evalue",0.001);
