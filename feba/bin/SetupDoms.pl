#!/usr/bin/perl -w
use Getopt::Long;
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for NewerThan()

{
    my $usage = "SetupDoms.pl [ -fast ] [ -hmmscan hmmscan ] nickname1 ... nicknameN\n";
    my $fast;
    my $hmmscan = "$Bin/hmmscan";
    GetOptions('hmmscan' => \$hmmscan, 'fast' => \$fast) || die $usage;
    my @orgs = @ARGV;
    die "Not executable: $hmmscan\n" unless -x $hmmscan;
    die "No such file: pfam.hmm\n" unless -e "pfam.hmm";
    die "No such file: tigrfam.hmm\n" unless -e "tigrfam.hmm";
    die $usage unless scalar(@orgs) > 0;

    my @to_build = ();
    foreach my $org (@orgs) {
	die "No g/$org directory" unless -d "g/$org";
        die "No aaseq file" unless -e "g/$org/aaseq";
        unless (NewerThan("g/$org/aaseq2","g/$org/aaseq")) {
            print STDERR "Making g/$org/aaseq2\n";
            system("./addOrgName.pl $org < g/$org/aaseq > g/$org/aaseq2") == 0 || die $!;
        }
        unless (NewerThan("g/$org/pfam.dom", "g/$org/aaseq")
            && NewerThan("g/$org/tigrfam.dom", "g/$org/aaseq")) {
            push @to_build, $org;
	    print STDERR "To build: g/$org/pfam.dom and g/$org/tigrfam.dom\n";
        }
    }

    if (scalar(@to_build) > 0) {
	if (defined $fast) {
	    print STDERR "Skipping HMMer for " . scalar(@to_build) . " directories\n";
	} else {
	    print STDERR "Running HMMer for " . scalar(@to_build) . " directories\n";
	    my $cmdsfile = "/tmp/SetupDoms.$$.cmds";
	    open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
	    foreach my $org (@to_build) {
		print CMDS "$hmmscan --cut_tc --domtblout g/$org/pfam.dom -o g/$org/pfam.hmmscan pfam.hmm g/$org/aaseq2 >& g/$org/pfam.log\n";
		print CMDS "$hmmscan --cut_tc --domtblout g/$org/tigrfam.dom -o g/$org/tigrfam.hmmscan tigrfam.hmm g/$org/aaseq2 >& g/$org/tigrfam.log\n";
	    }
	    close(CMDS) || die "Error writing to $cmdsfile";
	    system("$Bin/submitter.pl $cmdsfile") == 0 || die "Jobs failed";
	}
    }

    foreach my $org (@orgs) {
	foreach my $dom (qw{pfam tigrfam}) {
	    if (!-e "g/$org/$dom.dom") {
		print STDERR "Skipping $dom for $org, no hmmer table\n";
	    } else {
                print STDERR "Reading g/$org/$dom.dom\n";
		open(IN, "<", "g/$org/$dom.dom");
		open(OUT, ">", "g/$org/$dom.tab");
		print OUT join("\t", "locusId", "domainId", "domainName", "begin", "end", "score", "evalue")."\n";
		while(<IN>) {
		    chomp;
		    next if m/^#/;
		    my @F = split / +/, $_;
		    my ($domainName, $domainId, undef, $locusId, undef, $qlen,
			$evalueAll, $scoreAll, $biasAll,
			undef, undef,
			$evalue, undef, $score, $bias, $hmmBeg, $hmmEnd, $begin, $end) = @F;
		    die "Invalid locusId $locusId in $_" unless $locusId =~ m/:/;
		    die "Invalid qlen $qlen in $_" unless $qlen =~ m/^\d+$/;
		    die "Invalid begin or end $begin $end in $_" unless $begin =~ m/^\d+$/ && $end =~ m/^\d+$/;
		    die "Invalid score $score in $_: $score" unless  $score =~ m/^[0-9]+[.]?\d*$/;
		    die "Invalid begin end qlen $begin $end $qlen in $_" unless $begin >= 1 && $end >= $begin && $end <= $qlen;
		    die "Cannot parse evalue $evalue in $_" unless $evalue =~ m/^\d+[.]?\d*e?-?\d*$/;
		    $locusId =~ s/^.*://; # remove organism prefix
		    print OUT join("\t", $locusId, $domainId, $domainName, $begin, $end, $score, $evalue)."\n";
		}
		close(IN) || die "Error reading g/$org/$dom.dom";
		close(OUT) || die "Error writing g/$org/$dom.tab";
		print STDERR "Wrote g/$org/$dom.tab\n";
	    }
	}
    }

    foreach my $org (@orgs) {
        print STDERR "$Bin/RunRapSearch.pl $org\n";
        if (defined $fast) {
            print STDERR "Skipping...\n";
        } else {
            system("$Bin/RunRapSearch.pl",$org) == 0 || die "Error for rap search on $org";
        }
        print STDERR "Running KeggBestHit.pl on $org\n";
        system("$Bin/KeggBestHit.pl $org > g/$org/besthit.kegg") == 0
            || die "Error for KeggBestHit.pl on $org";
        unless(NewerThan("g/$org/seedanno.tab", "g/$org/aaseq")) {
          print STDERR "Running seed on $org\n";
          system("$Bin/run_seed.pl", $org);
        }
    }
}
