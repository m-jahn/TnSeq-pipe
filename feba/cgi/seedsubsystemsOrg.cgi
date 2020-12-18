#!/usr/bin/perl -w

#######################################################
## seedsubsystemsOrg.cgi -- SEED subsystems in an organism
##
## Copyright (c) 2017 University of California
##
## Author: Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

print
  header,
  Utils::start_page("SEED Subsystems in $orginfo->{$orgId}{genome}"),
  h2("SEED Subsystems in",
     a({ -href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome})),
  p("Subsystems are shown if they have at least one gene assigned to them. Because many roles are assigned to more than one subsystem, the subsystem may not be present even if the gene annotation(s) are correct.");

my $subsys = $dbh->selectall_arrayref(qq{ SELECT toplevel, category, subsystem, COUNT(DISTINCT locusId) AS nGenes
                                          FROM SEEDRoles JOIN SEEDAnnotationToRoles USING (seedrole)
                                          JOIN SeedAnnotation USING (seed_desc)
                                          WHERE orgId = ?
                                          GROUP BY subsystem
                                          ORDER BY toplevel, category, subsystem },
                                      { Slice => {} }, $orgId);

my @th = map th($_), ( "Category", "Subsystem", "#Genes" );
my @trows = ( \@th );
my $lastCat = undef;
foreach my $row (@$subsys) {
  my $subsys = $row->{subsystem};
  my $nice = $subsys; $nice =~ s/_/ /g;
  my $catShow = $row->{category} eq "" ? $row->{toplevel} : "$row->{toplevel} : $row->{category}";
  my @trow = ( td( defined $lastCat && $catShow eq $lastCat ? "&nbsp;" : $catShow ),
               td( a( { -href => "seedsubsystem.cgi?orgId=$orgId&subsystem=$subsys" }, $nice) ),
               td( { -style => "text-align: right;"}, $row->{nGenes} ) );
  $lastCat = $catShow;
  push @trows, \@trow;
}
  @trows = map Tr( { -valign => "top"}, @$_ ), @trows;
print
  table({ cellspacing => 0, cellpadding => 3}, @trows);

Utils::endHtml($cgi);
exit(0);

