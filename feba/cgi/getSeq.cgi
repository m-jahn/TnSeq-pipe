#!/usr/bin/perl -w
#######################################################
## getSeq.cgi
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Julie Jeong (jj326@berkeley.edu) and Morgan Price
#######################################################
use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use Bio::SeqIO;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();



# read the input information

my $orgId = $cgi->param('orgId') || "";
my $locusId = $cgi->param('locusId') || die "no locusId found\n";

my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $filename = $procId . $timestamp;
my $tmpDir = Utils::tmp_dir();
my $seqFile = "$tmpDir/$filename.fasta";
my $blast = '../bin/blast/blastall';
my $myDB = Utils::blast_db(); # sequence database
my $blastOut = "$tmpDir/$filename.blast.out";
my $blastSort = "$tmpDir/$filename.blast.sort";
my $seq;
my $dbh = Utils::get_dbh(); # gene database

print $cgi->header;

Utils::fail($cgi,"No orgId or no locusId specified") unless $orgId ne "" && $locusId ne "";

# extract sequence for the given gene
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
Utils::fail($cgi,"unknown gene") unless defined $gene->{locusId};
Utils::fail($cgi,"sequence information is only available for protein-coding genes.") unless $gene->{type} == 1;

my $fastacmd = '../bin/blast/fastacmd';
my $id = join(":",$orgId,$locusId);
system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!"; #error

my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
$seq = $in->next_seq()->seq;
$seq =~ s/(.{60})/$1\n/gs;

unlink($seqFile) || die "Error deleting $seqFile: $!";

my $showId = $gene->{sysName} || $gene->{locusId};
my $orginfo = Utils::orginfo($dbh);
my $title = "Sequence of $showId in $orginfo->{$orgId}{genome}",;
print
    start_html( -title => $title,
		-style => {-code => $style},
		-author=>'jj326ATberkeley.edu',
		-meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),
    h2("Sequence of $showId in " . $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),),
    pre(">$showId $gene->{desc} ($orginfo->{$orgId}{genome})\n$seq"),
    p(a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, "Show fitness")),
    p(a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Check homologs"));


Utils::endHtml($cgi);

exit 0;
