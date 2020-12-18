#!/usr/bin/perl -w

# Given the results of all-vs-all BLASTp on a database with sequence identifiers of the form
# organism:locus, make a table of "orthologs"

# Throughout, "locusId" is of the form organism:locus, not the second half only

use strict;
use Getopt::Long;
sub ReadBlastp($$);
sub BBH($$); # find bidirectional best BLAST hits with 80% coverage both ways

my $coverage = 0.8; # default requirement for coverage both ways
my $usage = "Usage: bbh.pl [-coverage $coverage]\n"
    . "        -blastp blastp_file -out out\n"
    . "   blastp_file must be in blast tab-delimited format, and locus identifiers must have the form\n"
    . "   taxonomy:locus\n"
    . "   Writes to out.bbh, with one line per gene, and out.scores, with one line per pair of BBHs\n";

my ($SUBJECT,$SCORE,$QBEGIN,$QEND,$SBEGIN,$SEND) = 0..5;

{
    my $compare = 0;
    my $blastpFile;
    my $prefix = undef;

    (GetOptions('coverage=f' => \$coverage,
		'blastp=s' => \$blastpFile,
		'out=s' => \$prefix)
     && @ARGV==0 && defined $blastpFile && defined $prefix) || die $usage;

    # taxa: taxonomyId => locusId => 1
    # hits: locusId => list of [subject,score,qbegin,qend,sbegin,send]
    my $taxa = {};
    my $hits = ReadBlastp($blastpFile,$taxa);
    print STDERR "Read hits for " . scalar(keys %$taxa) . " taxa\n";
    my $bbh = BBH($taxa,$hits);

    my @taxa = sort {$a cmp $b} (keys %$taxa);
    open(OUT, ">", "$prefix.bbh") || die "Cannot write to $prefix.bbh";
    open(SCORES, ">", "$prefix.scores") || die "Cannot write to $prefix.scores";
    print OUT join("\t","locusId","taxonomyId","")
	. join("\t", map { "orth" . $_ . ".locusId" } @taxa)."\n";
    print SCORES join("\t","tax1","locus1","tax2","locus2","ratio")."\n";

    foreach my $qtax (@taxa) {
	foreach my $query (keys %{ $taxa->{$qtax} }) {
	    my $q_bbh = $bbh->{$query};
	    my @orths = ();
	    @orths = map { exists $q_bbh->{$_} ? $q_bbh->{$_} : "" } @taxa;
	    print OUT join("\t",$query,$qtax,@orths)."\n";
	    my %maxscore = (); # query => subject => max score
	    foreach my $row (@{ $hits->{$query} }) {
		my ($subject,$score) = @$row;
		$maxscore{$subject} = $score if !exists $maxscore{$subject} || $maxscore{$subject} < $score;
	    }
	    if (!exists $maxscore{$query} && scalar(keys %$q_bbh) > 0) {
		print STDERR "Warning: BBH but no self hit for $query in $qtax\n";
	    } else {
		while (my ($stax,$subject) = each %$q_bbh) {
		    print SCORES join("\t", $qtax, $query, $stax, $subject, $maxscore{$subject} / $maxscore{$query})."\n";
		}
	    }
	}
    }
    close(OUT) || die "Error writing to $prefix.bbh";
    close(SCORES) || die "Error writing to $prefix.scores";
}

# could get memory intensive as #genomes increases
sub ReadBlastp($$) {
    my ($file,$taxa,$orth) = @_;
    my %hits = (); # locusId => list of [subject,score,qBeg,qEnd,sBeg,sEnd]
    open(BLASTP, "<", $file) || die "Cannot read $file";

    while (my $line = <BLASTP>) {
	chomp $line;
	my ($queryspec,$subjectspec,$identity,$aln_length,$mismatch,$gaps,$qBeg,$qEnd,$sBeg,$sEnd,$evalue,$score) = split /\t/, $line;
	die "Wrong number of columns in\n$line\n..." unless defined $score && $score =~ m/^ *[0-9.e+]+$/;
	$score =~ s/^ +//;
	my ($qtax,$query2) = split /:/, $queryspec; # to get taxId and locusId
	my ($stax,$subject2) = split /:/, $subjectspec;
	die "Invalid query $queryspec" unless defined $query2;
	die "Invalid subject $subjectspec" unless defined $subject2;
	my $query = "$qtax:$query2";
	my $subject = "$stax:$subject2";
	$taxa->{$qtax}{$query} = 1;
	$taxa->{$stax}{$subject} = 1;
	die $line unless defined $sEnd && $sEnd =~ m/^\d+$/;
	push @{ $hits{$query} }, [$subject,$score,$qBeg,$qEnd,$sBeg,$sEnd];
    }

    close(BLASTP) || die "Error reading $file";
    return(\%hits);
}

sub BBH($$) {
    my ($taxa,$hits) = @_;
    my %besthits = (); # query -> tax -> hit; does not include the coverage requirement
    my %len = (); # query -> length, based on self hit
    my %bbh = (); # locusId => taxonomyId => orthologId if it exists

    # reverse taxa to get the map of locusId => taxonomyId
    my %locusTax = ();
    while (my ($tax,$hash) = each %$taxa) {
	foreach my $locus (keys %$hash) {
	    $locusTax{$locus} = $tax;
	}
    }

    while (my ($query,$hitlist) = each %$hits) {
	my $qtax = $locusTax{$query};
	foreach my $row (@$hitlist) {
	    my ($subject,$score,$qBeg,$qEnd,$sBeg,$sEnd) = @$row;
	    if ($query eq $subject) {
		$len{$query} = $qEnd if !exists $len{$query} || $len{$query} < $qEnd;
	    } else {
		my $stax = $locusTax{$subject};
		next if $qtax eq $stax;
		$besthits{$query}{$stax} = $row unless
		    exists $besthits{$query}{$stax} && $besthits{$query}{$stax}[$SCORE] >= $score;
	    }
	}
    }
    # and then combine besthits and check coverage
    while (my ($query,$bhits) = each %besthits) {
	my $qtax = $locusTax{$query};
	while (my ($stax,$row) = each %$bhits) {
	    my ($subject,$score,$qBeg,$qEnd,$sBeg,$sEnd) = @$row;
	    die unless defined $qEnd;
	    if (!defined $len{$query}) {
		print STDERR "No self-hit for $query\n";
		next;
	    }
	    if (!defined $len{$subject}) {
		print STDERR "No self-hit for $subject\n";
		next;
	    }
	    if ($qEnd-$qBeg+1 >= $coverage * $len{$query}
		&& $sEnd-$sBeg+1 >= $coverage * $len{$subject}
		&& exists $besthits{$subject}{$qtax}
		&& $besthits{$subject}{$qtax}[$SUBJECT] eq $query) {
		$bbh{$query}{$stax} = $subject;
		$bbh{$subject}{$qtax} = $query;
	    }
	}
    }
    return \%bbh;
}
