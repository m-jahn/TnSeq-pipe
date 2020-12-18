#!/usr/bin/perl -w

#######################################################
## downloadReanno.cgi
##
## Copyright (c) 2018 University of California
##
## Authors:
## Morgan Price
#######################################################

# Download a table of reannotated proteins and their sequences, either
# for one organism or all of them.
#
# Optional CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use Time::HiRes qw(gettimeofday);

use lib "../lib";
use Utils;
use FEBA_Utils; # for ReadFasta()

my $cgi = CGI->new;
my $orgId = $cgi->param('orgId') || "";
my $dbh = Utils::get_dbh();
my $myDB = Utils::blast_db();

my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId $orgId" if $orgId ne "" && !exists $orginfo->{$orgId};

# Ignore orgId if empty
my $cmpfunc = $orgId eq "" ? "<>" : "=";
my $reanno = $dbh->selectall_arrayref("SELECT * from Reannotation JOIN Gene USING (orgId,locusId) WHERE orgId $cmpfunc ?",
                                      { Slice => {} }, $orgId);
my $filename = $orgId ? "reanno_$orgId.tsv" : "reanno.tsv";

# Fetch the sequences
my $aaseqs = {}; # lcl|orgId:locusId => sequence
if (@$reanno > 0) {
  my $procId = $$;
  my $timestamp = int (gettimeofday * 1000);
  my $pre = Utils::tmp_dir() . "/" . $procId . $timestamp;
  my $listFile = "$pre.list";
  my $faaFile = "$pre.faa";
  open(LIST, ">", $listFile) || die "Cannot write to $listFile";
  foreach my $row (@$reanno) {
    print LIST join(":", $row->{orgId}, $row->{locusId}) . "\n";
  }
  close(LIST) || die "Error writing to $listFile";
  my $fastacmd = '../bin/blast/fastacmd';
  system($fastacmd, '-d', $myDB, '-i', $listFile, '-o', $faaFile) == 0
    || die "Error running $fastacmd -- $!";
  $aaseqs = FEBA_Utils::ReadFasta($faaFile);
  unlink($listFile);
  unlink($faaFile);
}

print "Content-Type:text/tab-separated-values\n";
print "Content-Disposition: attachment; filename=$filename\n\n";

print join("\t", qw/organism orgId locusId sysName gene reannotation comment original_desc aaseq/)."\n";
foreach my $row (@$reanno) {
  my $id = "lcl|" . $row->{orgId} . ":" . $row->{locusId};
  die "No sequence for $id" unless exists $aaseqs->{$id};
  print join("\t", $orginfo->{ $row->{orgId} }{genome}, $row->{orgId},
             $row->{locusId}, $row->{sysName} || "", $row->{gene} || "",
             $row->{new_annotation}, $row->{comment}, $row->{desc},
             $aaseqs->{$id})."\n";
}
$dbh->disconnect();
