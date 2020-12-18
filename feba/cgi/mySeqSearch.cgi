#!/usr/bin/perl -w
#######################################################
## mySeqSearch.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################
#
# Required CGI parameters:
# Either query (the sequence, a.a. by default)
# or orgId and locusId (if coming from a page for that gene)
#
# Optional CGI parameters in query mode:
# qtype -- protein or nucleotide (default is protein)
#
# Optional CGI parameters in either mode:
# numHit -- how many hits to show

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use Bio::SeqIO;

use lib "../lib";
use Utils;
use IO::Handle;
use URI::Escape;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $nCPU = $ENV{MC_CORES} || 8;

# read the input information
my $query = $cgi->param('query') || "";
my $qtype = $cgi->param('qtype') || "protein";
my $numHit = $cgi->param('numHit') || 100;
my $locusSpec = $cgi->param('locusId') || "";
my $orgId = $cgi->param('orgId') || "";
my $dbh = Utils::get_dbh();

# set up $seq and $seqFile
my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $filename = $procId . $timestamp;
my $tmpDir = Utils::tmp_dir();
my $seqFile = "$tmpDir/$filename.fasta";
my $blast = '../bin/blast/blastall';
my $myDB = Utils::blast_db();
my $blastOut = "$tmpDir/$filename.blast.out";
my $blastSort = "$tmpDir/$filename.blast.sort";
my $seq;
my $locusShow;
my $fastacmd = '../bin/blast/fastacmd';
die "No such executable: $fastacmd" unless -x $fastacmd;

print $cgi->header; # must be printed before using fail()

my $def = "";
my $gene;
if ($query =~ m/[A-Za-z]/) {
    # parse and write the input sequence
    $seq = "";
    my @lines = split /[\r\n]+/, $query;
    $def = shift @lines if @lines > 0 && $lines[0] =~ m/^>/;
    $def =~ s/^>//;
    foreach (@lines) {
        s/[ \t]//g;
        s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
        next if $_ eq "//";
        Utils::fail($cgi,"Error: more than one sequence was entered.") if m/^>/;
        Utils::fail($cgi,"Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
        s/[*]/X/g;
        $seq .= uc($_);
    }
    $def = substr($seq,0,20) . "..." if $def eq "";

    my $id = "query";
    open(FAA,">",$seqFile) || die "Cannot write fasta file";
    print FAA Utils::formatFASTA($id,$seq);
    close(FAA) || die "Error writing fasta file";    
} elsif ($locusSpec ne "") {
    Utils::fail($cgi,qq($locusSpec is invalid. Please enter correct gene name!)) unless ($locusSpec =~ m/^[A-Za-z0-9_.-]*$/);

    # extract sequence for the given gene
    $gene = $dbh->selectrow_hashref("SELECT * from Gene where orgId=? AND locusId=?", {}, $orgId, $locusSpec);

    Utils::fail($cgi, "no such gene: $locusSpec in $orgId") unless defined $gene->{locusId};
    Utils::fail($cgi, "homology search is only available for protein-coding genes") unless $gene->{type} == 1;
    my $id = join(":",$orgId,$locusSpec);
    $locusShow = $gene->{gene} || $gene->{sysName} || $gene->{locusId};
    $def = $locusShow;
    system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!";
    my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
    $seq = $in->next_seq()->seq;
} else {
    Utils::fail($cgi, "No sequence or gene specified!");
}

# print start of page
my $orginfo = Utils::orginfo($dbh);
my $title = "";
my $tabs = "";
my $orth;
my $showOrth = $orgId && $locusSpec;
if ($showOrth) {
    $title = "Homologs of $locusShow from $orginfo->{$orgId}{genome}";
    $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusSpec,0,1,"homo");
    $orth = $dbh->selectall_hashref("SELECT * FROM Ortholog WHERE orgId1 = ? AND locusId1 = ?",
                             'orgId2',{Slice=>{}}, $orgId, $locusSpec);
} else {
    my $qlen = length($seq);
    my $qchar = $qtype eq "protein" ? "a.a." : "nt.";
    $title = "Blast Results for " . substr($seq,0,20) . "... ($qlen $qchar)";
    $tabs = '<div id = "ntcontent">';
}
my $start = Utils::start_page($title);

print $start, $tabs, $cgi->h2($title);
print p(Utils::gene_link($dbh, $gene, "lines")) if $showOrth;
print "\n";
autoflush STDOUT 1; # so header shows while blast is being run

# run blast
if ($qtype eq "nucleotide") {
    Utils::fail($cgi,qq($query is invalid. Please enter nucleotide sequence or choose sequence type as protein!)) unless ($seq =~ m/^[ATCGatcg]*$/);
    system($blast,'-p','blastx','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8','-a',$nCPU)==0 || die "Error running blastx, are you sure this is a nucleotide query: $!";
} elsif ($qtype eq "protein") {
    Utils::fail($cgi,qq($query is invalid. Please enter correct protein sequence!)) unless ($seq =~ m/^[A-Za-z]*$/);
    system($blast,'-p','blastp','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8','-a',$nCPU)==0 || die "Error running blastp, are you sure this is a protein query: $!";
} else {
    die "Unknown query type $qtype";
}

# parse and report the blast result:
# blast output fields: (1)queryId, (2)subjectId, (3)percIdentity, (4)alnLength, (5)mismatchCount, (6)gapOpenCount, (7)queryStart, (8)queryEnd, (9)subjectStart, (10)subjectEnd, (11)eVal, (12)bitScore
# sort the blast result by bit score, E-value, and percent identity
system('sort','-k1,1','-k12,12gr','-k11,11g','-k3,3gr',$blastOut,'-o',$blastSort)==0 || die "Error running sort: $!";

# output blast result

my $defencoded = uri_escape($def); # the defline may have wierd characters
my $aln_preURL = $query ? "showAlign.cgi?query=$defencoded&querySequence=$seq"
  : "showAlign.cgi?query=$orgId:$locusSpec";

open(RES,$blastSort) || die "Error reading $blastSort";
my @hits = (); # each row is a simple list of the fields in the blast ouptut
while(<RES>) {
  chomp;
  my @F = split /\t/, $_;
  die "Truncated blast output" unless @F == 12;
  push @hits, \@F;
}
close(RES) || die "Error reading $blastOut";
my $trunc = 0;
if (defined $numHit && $numHit > 0 && @hits > $numHit) {
  $trunc = 1;
  $#hits = $numHit - 1
}

if (@hits > 0) {
  print $cgi->p(($trunc ? "Top " : "") . scalar(@hits) . " hits from BLASTp (E < 0.01)");
  my @header = (a({-title => "Potential ortholog from bidirectional best hit?", -style => "color: black;"}, 'Orth'),
                'Species', 'Gene', 'Name', 'Description', 'Fitness',
                a({-title => "Percent identity", -style=>"color: black;"}, '%Id'),
                a({-title => "Percent coverage of query", -style=>"color: black;"}, 'Cov'));
  @header = map small($_), @header;
  # Hard code the widths to avoid annoying table resize during incremental rendering
  my @widths = ("4%", "23%", "14%", "5%", "37%", "6%","5%","5%");
  my @th = map th({width => $widths[$_]}, $header[$_]), 0..$#header;
  shift @th unless $showOrth;
  # hard-coding the table width prevents annoying resize during incremental updates
  # and the strong bottom border on the table is distracting during incremental updates so remove it
  print qq{<TABLE style="border-bottom: none;" cellspacing=0 cellpadding=3 width=99% >},
    Tr(@th),
    "\n";

  foreach my $hit (@hits) {
    my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,
        $queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitScore) = @$hit;
    $bitScore =~ s/ +//;
    my ($orgId,$locusId) = split /:/, $subjectId;
    my $cov = int(0.5 + 100*abs($queryEnd - $queryStart + 1)/length($seq));
    $percIdentity = int(0.5 + $percIdentity);
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                       undef, $orgId, $locusId);
    if (!defined $gene) {
      print "Warning! Unknown hit $orgId:$locusId<BR>";
      next;
    }
    my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $orgId, $locusId);
    my $showId = $gene->{sysName} || $gene->{locusId};
    my $seqlen = length($seq);
    my $aln_URL = "$aln_preURL&subject=$subjectId"
      unless $qtype eq "nucleotide";

    # Compute length of the subject
    open(my $fh, "-|", $fastacmd, '-d',$myDB, '-s', $subjectId)
      || die "Cannot run fastacmd";
    my $in = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    my $sseq = $in->next_seq()->seq;
    die "No sequence for subject $subjectId in $myDB" unless $sseq;
    my $slen = length($sseq);
    close($fh) || die "Error running $fastacmd";

    my $cov_title = "amino acids $queryStart:$queryEnd / $seqlen of query"
      . " are similar to a.a. $subjectStart:$subjectEnd / $slen of $showId";
    my $covShow = defined $aln_URL ? a({ title => $cov_title, href => $aln_URL }, $cov)
      : a({ title => $cov_title }, $cov);
    my @row = ($cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}->{genome}),
               Utils::gene_link($dbh, $gene, "name", "geneOverview.cgi"),
               $gene->{gene},
               Utils::gene_link($dbh, $gene, "desc", "domains.cgi"),
               $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId", title => $fittitle }, $fitstring ),
               $cgi->a({title=>"evalue: $eVal ($bitScore bits)"},$percIdentity),
               $covShow);
    my @td = map td({width => $widths[$_+1]}, $row[$_]), (0..$#row);
    if ($showOrth) {
      my $o = exists $orth->{$orgId} && $orth->{$orgId}{locusId2} eq $locusId ?
        a({ -title => "Potential ortholog (bidirectional best hit with high coverage)"}, "o")
          : '&nbsp;';
      unshift @td, td({width => $widths[0], style=>"text-align: center;"}, $o);
    }
    print $cgi->Tr({ -align => 'left', -valign => 'top', bgcolor=>'white' }, @td);
    print "\n";
  }
  print "</TABLE>";
} else {
  print $cgi->p("No hits found using BLASTp (E < 0.01)");
}

print qq[<br>Or search for homologs using
	<A HREF="http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>$def%0A$seq"
           TITLE="Find papers about this protein or its homologs">PaperBLAST</A>,
       <A HREF="https://iseq.lbl.gov/genomes/seqsearch?sequence=>$def%0A$seq">ENIGMA genome browser</A>,
       or
       <A HREF="http://www.microbesonline.org/cgi-bin/seqsearch.cgi?qtype=protein&query=$seq">MicrobesOnline</A>
       <BR>
       <BR>] unless $qtype eq "nucleotide";
$dbh->disconnect();
unlink($seqFile) || die "Error deleting $seqFile: $!";
unlink($blastOut) || die "Error deleting $blastOut: $!";
unlink($blastSort) || die "Error deleting $blastSort: $!";

Utils::endHtml($cgi);
exit 0;
