#!/usr/bin/perl -w
#######################################################
## exps.cgi -- search for experiments, or list all in one organism
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price and Victoria Lo
#######################################################
#
# Key parameters: orgId and query, for which organism and to look up experiments
#	At least one must be meaningful (present and not empty)
# OR, specify expGroup AND condition1. In this set up, condition1 may be empty,
#	but it must still be specified (and is used to restrict the results).
# OR, specify orgId AND expGroup without condition1.

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;
use HTML::Entities;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $expSpec = $cgi->param('query');
$expSpec = "" if !defined $expSpec;
$expSpec =~ s/ +$//;
$expSpec =~ s/^ +$//;
$expSpec =~s/[\"\n\r\;\\]//g; # ' is allowed, handled in matching_exps

my $expGroup = $cgi->param('expGroup');
my $condition1 = $cgi->param('condition1');

my $dbh = Utils::get_dbh();

# Redirect to org.cgi if displaying all experiments from one organism,
# or to orgAll.cgi if displaying all experiments
if ($orgId ne "" && !defined $expGroup && ($cgi->param("All experiments") || $expSpec eq "")) {
    print redirect(-url=>"org.cgi?orgId=$orgId");
    exit(0);
} elsif ($orgId eq "" && !$expSpec && !$expGroup && !$condition1) {
  print redirect(-url=>"orgAll.cgi");
  exit(0);
}

my $specShow = $expSpec;
if ($specShow eq "") {
    $specShow = defined $condition1 ? join(" ", encode_entities($expGroup), encode_entities($condition1))
      : encode_entities($expGroup);
    $specShow = "" if !defined $specShow;
}
my $start = Utils::start_page("Experiments for $specShow");
$expSpec = "" if $cgi->param("All experiments");

print $cgi->header, $start, '<div id="ntcontent">';

my $exps;

if (defined $expGroup && !defined $condition1){
    $exps = $dbh->selectall_arrayref(qq{SELECT * from Experiment WHERE expGroup = ? AND orgId = ? GROUP BY expGroup, condition_1, expName},
            { Slice => {} },
            $expGroup, $orgId);
} elsif (defined $expGroup && defined $condition1) {
    $exps = $dbh->selectall_arrayref(qq{SELECT * from Experiment WHERE expGroup = ? AND condition_1 = ? GROUP BY expGroup, condition_1, expName},
				    { Slice => {} },
				    $expGroup, $condition1);
    Utils::fail($cgi, "No experiments for specified group and condition") if scalar(@$exps) == 0;
} elsif ($orgId eq "" && defined $expSpec) {
    $exps = Utils::matching_exps($dbh, "", $expSpec);
} else {
    $exps = Utils::matching_exps($dbh, $orgId, $expSpec);
}

# Make sure both parameters are safe
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};

# sort experiments by organism, or by lc(condition)
if ($orgId eq "") {
    my @exps = sort { $orginfo->{$a->{orgId}}{genome} cmp $orginfo->{$b->{orgId}}{genome} } @$exps;
    $exps = \@exps;
}
if ($orgId && defined $expGroup) {
    my @exps = sort { lc($a->{condition_1}) cmp lc($b->{condition_1}) } @$exps;
    $exps = \@exps;
}

if (@$exps == 0) {
   print $cgi->h3(qq{No experiment found matching "$expSpec"});
} else {
  my $heading = "Experiments";
  $heading .= " in ". $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}") if $orgId ne "";
  if ($expSpec ne "") {
      $heading .= qq{ matching "$expSpec"};
  } elsif ($specShow ne "") {
      $heading .= " for $specShow";
  }
  print $cgi->h2($heading);

  my $exp1 = $exps->[0];
  if ($exp1->{expGroup} ne ""
      && ($expSpec ne "" || (defined $expGroup && defined $condition1))) {
      print p(a( { -href => "orthCond.cgi?expGroup=" . uri_escape($exp1->{expGroup})
                       . "&condition1=" . uri_escape($exp1->{condition_1}) },
                 "Or see specific phenotypes for $exp1->{expGroup} $exp1->{condition_1} across organisms"));
  }
}

  my @trows = ();


# remove organism column if already specified
if (defined $expGroup && defined $orgId && !defined $condition1){
   push @trows, $cgi->Tr({-valign => "top"}, $cgi->th(['Name', 'Group', 'Condition', 'Description' ]));
  foreach my $row (@$exps) {
    my @cond = map $row->{"condition_".$_}, (1..4);
    @cond = grep $_ ne "", @cond;
    push @trows, $cgi->Tr({-valign => "top"},
                          $cgi->td([$cgi->a({href => "exp.cgi?orgId=$row->{orgId}&expName=$row->{expName}"}, $row->{expName}),
                          $row->{expGroup},
                          join("; ", @cond),
                          $row->{expDesc} ]));
    }
} else {
  push @trows, $cgi->Tr({-valign => "top"}, $cgi->th([ 'Organism', 'Name', 'Group', 'Condition', 'Description' ]));
  foreach my $row (@$exps) {
      push @trows, $cgi->Tr({-valign => "top"},
             $cgi->td([ $cgi->a({href => "org.cgi?orgId=". $orginfo->{$row->{orgId}}->{orgId}}, "$orginfo->{$row->{orgId}}{genome}"),
	              $cgi->a({href => "exp.cgi?orgId=$row->{orgId}&expName=$row->{expName}"}, $row->{expName}),
		      $row->{expGroup}, $row->{condition_1}, $row->{expDesc} ]));
  }
}

print table({cellspacing => 0, cellpadding => 3}, @trows),
    "<br><br>";
print p("Or see all specific phenotypes in ",
        a({href => "spec.cgi?orgId=$orgId&expGroup=" . uri_escape($expGroup)}, $expGroup),
        " experiments for $orginfo->{$orgId}{genome}")
    if defined $expGroup && @$exps > 0 && $orgId ne "";

$dbh->disconnect();
Utils::endHtml($cgi);
