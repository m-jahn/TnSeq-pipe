#!/usr/bin/perl -w

#######################################################
## seedsubsystem.cgi -- overview of a SEED subsystem
##
## Copyright (c) 2017 University of California
##
## Author: Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# subsystem -- which subsystem
#
# Optional CGI parametrs:
# expName -- which experiment(s) to show
# expQuery -- which experiment(s) to add

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $subsystem = $cgi->param('subsystem') || die "No subsystem parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

my $nice = $subsystem; $nice =~ s/_/ /g;
print
  header,
  Utils::start_page("SEED Subsystem: $nice in $orginfo->{$orgId}{genome}"),
  h2("SEED Subsystem: $nice in",
     a({ -href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome}));

my $roles = $dbh->selectcol_arrayref("SELECT seedrole FROM SEEDRoles WHERE subsystem = ?",
                                     {}, $subsystem);

my @exps = (); # experiment objects in order
my @expNames = param('expName');

my %expNames = (); # expName to experiment object
foreach my $expName (@expNames) {
  next if $expName eq "";
  my $exp = $dbh->selectrow_hashref("SELECT * from Experiment WHERE orgId = ? AND expName = ?",
                                    {}, $orgId, $expName);
  die "No such experiment for $orgId: $expName" unless $exp;
  push @exps, $exp unless exists $expNames{$expName};
  $expNames{$expName} = $exp;
}

my $expQuery = param('expQuery');
if ($expQuery) {
  my @exps2 = @{ Utils::matching_exps($dbh,$orgId,$expQuery) };
  print p(qq{Sorry, no experiments match "$expQuery"})
    if @exps2 == 0;
  foreach my $exp (@exps2) {
    push @exps, $exp unless exists $expNames{$exp->{expName}};
    $expNames{$exp->{expName}} = $exp;
  }
}

if (@$roles == 0) {
  print p("No roles -- is this a valid subsystem?");
} else {
  # Experiment selector UI
  my @hidden = (hidden('orgId'), hidden('subsystem'));
  foreach my $exp (@exps) {
    push @hidden, hidden( -name => 'expName', -value => $exp->{expName}, -override => 1 );
  }
  print qq[<div style="position: relative;"><div class="floatbox">],
    p(start_form(-name => 'input', -method => 'GET', -action => 'seedsubsystem.cgi'),
      join("", @hidden),
      "Experiments:",
      textfield( -name => 'expQuery', -default => '', -override => 1, -size => 15, -maxLength => 100),
      "<BUTTON type='submit'>Add</BUTTON>",
      end_form),
    qq[</div></div>];
  # Columns are: role, gene, gene name, experiment column(s) if any, checkbox column
  my @th = map th($_), qw{Role Gene Name};
  foreach my $exp (@exps) {
    push @th, th(a({ -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}", -title => $exp->{expName} },
                   $exp->{expDesc}));
  }
  push @th, th("&nbsp;");
  my @trows = ( \@th );
  if (@exps > 0) {
    my $rmmark = "&#10799;"; # Unicode Character 'VECTOR OR CROSS PRODUCT' (U+2A2F)
    # Add remove-experiment controls
    my @trow = map td($_), ("","","");
    foreach my $exp (@exps) {
      my @args = ("orgId=$orgId", "subsystem=$subsystem");
      foreach my $exp2 (@exps) {
        push @args, "expName=$exp2->{expName}" unless $exp2->{expName} eq $exp->{expName};
      }
      my $rmURL = "seedsubsystem.cgi?" . join("&", @args);
      push @trow, td({ -align => "center" },
                     a({ -href => $rmURL, -title => "remove $exp->{expName}" }, $rmmark));
    }
    push @trow, td("");
    push @trows, \@trow;
  }
  foreach my $role (@$roles) {
    my $genes = $dbh->selectall_arrayref(qq{ SELECT * FROM SEEDAnnotationToRoles
                                             JOIN SEEDAnnotation USING (seed_desc)
                                             JOIN Gene USING (orgId, locusId)
                                             WHERE orgId = ? AND seedrole = ? },
                                         { Slice => {} }, $orgId, $role);
    my $showRole = $role;
    # Link to EC number hits from other sources
    $showRole =~ s!EC (\d+[.]\d+[.]\d+[.]\d+)!<A HREF="myFitShow.cgi?orgId=$orgId&gene=ec:$1">EC $1</A>!g;
    if (@$genes == 0) {
      my @trow = map td($_), ( $showRole, "No genes", "");
      foreach my $exp (@exps) { push @trow, td(""); }
      push @trow, td("");
      push @trows, \@trow;
    } else {
      my $first = 1;
      foreach my $gene (@$genes) {
        my @trow = map td($_), ( $first ? $showRole : "&nbsp;",
                                 Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi"),
                                 $gene->{gene} || "");
        foreach my $exp (@exps) {
          my $showfit = "&nbsp;";
          my ($fit,$t) = $dbh->selectrow_array("SELECT fit,t FROM GeneFitness WHERE orgId = ? AND locusId = ? AND expName = ?",
                                               {}, $orgId, $gene->{locusId}, $exp->{expName});
          if (defined $fit) {
            my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$gene->{locusId}&expName=$exp->{expName}";
            $showfit = a({ -href => $strainUrl,
                           -title => "t = " . sprintf("%.1f", $t),
                           -style => "color:rgb(0,0,0)" },
                         sprintf("%.1f", $fit));
          }
          push @trow, td({ -bgcolor => Utils::fitcolor($fit) }, $showfit);
        }
        push @trow, td(checkbox('locusId', 1, $gene->{locusId}, ''));
        push @trows, \@trow;
        $first = 0;
      }
    }
  }
  @trows = map Tr( { -valign => "top", -align => "left" }, @$_ ), @trows;
  print start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
    hidden('orgId', $orgId),
      table({ cellspacing => 0, style => "text-align: center;", cellpadding => 3}, @trows),
	p(submit(-class=>"heatmap", -name=>"Heatmap of selected genes")),
          end_form;
}

print
  p( start_form(-name => 'orgselect', -method => 'GET', -action => 'seedsubsystem.cgi'),
     hidden( -name => 'subsystem', -value => $subsystem, -override => 1),
     "Or select organism:",
     Utils::OrgSelector($orgId, $orginfo),
     "<button type='submit'>Go</button>",
     end_form ),
  p(a( { -href => "http://pubseed.theseed.org/FIG/seedviewer.cgi?page=Subsystems&subsystem=$subsystem" },
       "Or view this subsystem at the SEED" ),
   "(slow)");

Utils::endHtml($cgi);
exit(0);

