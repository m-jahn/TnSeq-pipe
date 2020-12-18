#!/usr/bin/perl -w

#######################################################
## keggmap.cgi
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Show a KEGG map
#
# CGI parameters -- mapId, such as map00010 or 00010
# optional:
#	orgId, to specify which enzymes to mask out and which to search for members in
#	expName, to specify which experiment(s) to color code by
#	or expquery, to search for experiments to use
#	ec -- a comma-delimited list of EC #s to highlight

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

# LoadKEGGMap(dbh, kegg id) => map data structure, including
# mapId
# mapdesc
# links -- a list of [coord, url] -- these are the non-enzyme links only (as they do not get shaded)
# ecnums -- a list of hash of ecnum, ecdesc, coord, url, and optionally
#	mask=1 (to grey out), or else colors = ref to list of colors
#	From this routine, mask is never set, and
#	the urls are either to kegg (ww.genome.jp) or they search for the EC #
sub LoadKEGGMap($$);

# map data structure => outputs HTML for the map image and overlays
#	uses the above "map" data structure
#	also uses an optional highlight attribute for each ecnum entry
sub DrawKEGGMap($);

my $cgi=CGI->new;
my $mapId = $cgi->param('mapId');
my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;
my $expquery = $cgi->param('expquery');
my $ecHighlight = $cgi->param('ec') || "";
my @ecHighlight = split /,/, $ecHighlight;
my %ecHighlight = map { $_ => 1 } @ecHighlight;

my @expNames = $cgi->param('expName');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown orgId $orgId" if $orgId && !exists $orginfo->{$orgId};

die "Must specify mapId" unless defined $mapId;
$mapId =~ s/^map//;
die "Invalid mapId" unless $mapId =~ m/^\d+$/;

my $map = &LoadKEGGMap($dbh, $mapId);

my $title = $map->{mapdesc};
my $title2 = $title;
if ($orgId) {
    $title .= " in " . $orginfo->{$orgId}{genome} ;
    $title2 .= " in " . a({href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome});
}
my $start = Utils::start_page($title);
print
    header,
    $start, '<div id="ntcontent">',
    h2($title2);

if ($orgId && defined $expquery) {
    my @exps = @{ Utils::matching_exps($dbh,$orgId,$expquery) };
    Utils::fail($cgi, qq{No experiment matching "$expquery"}) if @exps == 0;
    if (@exps == 1) {
        @expNames = map $_->{expName}, @exps;
    } else {
        # choose from among the matches
        my @trows = ();
        my @headings = qw{&nbsp; name group condition description};
        push @trows, $cgi->Tr({-valign => 'top', -align => 'center'}, $cgi->th(\@headings));
        foreach my $exp (@exps) {
            push @trows, $cgi->Tr(
                {-valign => 'top', -align => 'left'},
                $cgi->td([ qq{<input type="checkbox" name="expName" value="$exp->{expName}"  >},
                           $cgi->a({href => "exp.cgi?orgId=$orgId&$exp->{expName}"}, $exp->{expName}),
                           $exp->{expGroup}, $exp->{condition_1}, $exp->{expDesc} ]));
        }
        print
            p("Select experiment(s) to highlight"),
            start_form(-name => 'input', -method => 'GET', -action => 'keggmap.cgi'),
            hidden( -name => 'mapId', -value => $mapId, -override => 1),
            hidden( -name => 'orgId', -value => $orgId, -override => 1),
            table( {cellpadding => 3, cellspacing => 0}, @trows),
            submit('Go'),
            end_form;
        Utils::endHtml($cgi); #exits
    }
}

my %ecs = map { $_->{ecnum} => 1 } @{ $map->{ecnums} };
my @ecs = keys(%ecs);

my @expShow = (); # HTML for each highlighted experiment
my $expspecParams = ""; # i.e., "&expName=set1H5&expName=set1H6"

if ($orgId) {
    my $expinfo = Utils::expinfo($dbh,$orgId);
    foreach my $expName (@expNames) {
        die "No such experiment: $expName" unless exists $expinfo->{$expName};
        push @expShow, a({href => "exp.cgi?orgId=$orgId&expName=$expName",
                          title => "$expName: $expinfo->{$expName}{expDescLong}"},
                         $expinfo->{$expName}{expDesc});
        $expspecParams .= "&expName=$expName";
    }
    my $expSpec = join(",", map { "'" . $_ . "'" } @expNames);

    # grey out absent EC#s, color code (if @expNames), and link to actual genes if present
    my $ecGenes = Utils::EcToGenes($dbh, $orgId, \@ecs);
    foreach my $row (@{ $map->{ecnums} }) {
        my $ecnum = $row->{ecnum};
        $row->{highlight} = 1 if exists $ecHighlight{$ecnum};
        if (exists $ecGenes->{$ecnum}) {
            my @locusIds = sort keys %{ $ecGenes->{$ecnum} };
            my @locusSpecs = map "locusId=$_", @locusIds;
            if (@locusSpecs > 1) {
                $row->{url} = "genesFit.cgi?orgId=$orgId&" . join("&",@locusSpecs);
            } else {
                $row->{url} = "singleFit.cgi?orgId=$orgId&" . $locusSpecs[0];
            }
            # and color code
            if (@expNames) {
                my $locusIn = join(",", map "'" . $_ . "'",  @locusIds);
                my $fitrows = $dbh->selectall_arrayref(
                    qq{ SELECT locusId, avg(fit) FROM GeneFitness
                        WHERE locusId IN ( $locusIn ) AND expName IN ( $expSpec ) AND orgId = ?
                        GROUP BY locusId
                        ORDER BY locusId },
                    {}, $orgId);
                my %fit = map { $_->[0] => $_->[1] } @$fitrows;
                # note undef => grey
                my @colors = map Utils::fitcolor($fit{$_}), @locusIds;
                $row->{colors} = \@colors;

                # and update hover text
                $row->{ecdesc} =~ s/[.]$//; # removing trailing period
                my @sorted = sort { $a <=> $b} values %fit;
                if (scalar(keys %fit) == 1) {
                    $row->{ecdesc} = sprintf("%s, fitness %.1f",
                                             $row->{ecdesc}, $sorted[0]);
                } elsif (scalar(keys %fit) > 1) {
                    $row->{ecdesc} = sprintf("%s, fitness %.1f to %.1f",
                                             $row->{ecdesc}, $sorted[0], $sorted[@sorted-1]);
                } else { # genes but no values
                    $row->{ecdesc} .= ", no fitness data";
                }
            }
        } else {
            $row->{mask} = 1;
            $row->{ecdesc} =~ s/[.]$//; # removing trailing period
            $row->{ecdesc} .= " (no genes)";
        }
    }
    # pathway links to remember orgId etc.
    foreach my $row (@{ $map->{links} }) {
        # entry 1 is the url
        $row->[1] .= "&orgId=$orgId".$expspecParams if $row->[1] =~ m/^keggmap.cgi/;
    }
}

print
    p(a( {href => "keggmaplist.cgi?orgId=$orgId".$expspecParams}, "Browse metabolic maps")),
    p( start_form(-name => 'orgselect', -method => 'GET', -action => 'keggmap.cgi'),
       hidden( -name => 'mapId', -value => $mapId, -override => 1),
       "Select organism:",
       Utils::OrgSelector($orgId, $orginfo),
       "<button type='submit'>Go</button>",
       end_form );

if (@expShow) {
    print p("Colored by",
            (@expShow > 0 ? "average" : ""),
            "fitness in:",
            join(", ", @expShow),
            "<BR/>",
            a({href => "keggmap.cgi?mapId=$mapId&orgId=$orgId" }, "(clear)"));
} elsif ($orgId) {
    print
        p( start_form('-name' => 'expsearch', -method => 'GET', -action => 'keggmap.cgi'),
           hidden( -name => 'mapId', -value => $mapId, -override => 1),
           hidden( -name => 'orgId', -value => $orgId, -override => 1),
           "Color by fitness in:",
           qq{ <input type="text" name="expquery" size="20" maxlength="100" /> },
           qq{ <button type='submit'>Go</button> },
           end_form );
}

&DrawKEGGMap($map);

print
    h3("Help"),
    ul(li("Enzyme classification numbers are greyed out if $orginfo->{$orgId}{genome} is not predicted to contain this reaction."),
       li("These predictions are often incorrect, especially if the EC number is incomplete (as with '2.7.1.-')."),
       li("Gene-enzyme associations are based on the last (2011) public release of the",
          a( {href => "http://www.genome.jp/kegg/"}, "Kyoto Encyclopedia of Genes and Genomes"),
          "(KEGG) and also on",
          a( {href => "http://www.jcvi.org/cgi-bin/tigrfams/index.cgi"}, "TIGRFAMS"),
          "and",
          a( {href => "http://theseed.org"}, "SEED") . ".") )
    if $orgId ne "";


$dbh->disconnect();
Utils::endHtml($cgi);

# returns a hashref that includes mapId, mapdesc, links, ecnums
sub LoadKEGGMap($$) {
    my ($dbh,$mapid) = @_;
    my ($mapdesc) = $dbh->selectrow_array("SELECT title FROM KEGGMap where mapId = ?",
                                      {}, $mapId);
    die "Unknown map id $mapId" unless defined $mapdesc;
    my $mapobjs = $dbh->selectall_arrayref("SELECT objectId,type,coord,url FROM KEGGConf WHERE mapId = ?",
                                           {}, $mapId);
    my @links = ();
    my @ecnums = ();
    foreach my $row (@$mapobjs) {
        my ($objectId,$type,$coord,$url) = @$row;
        $url = "http://www.genome.jp$url"; # default is, link to KEGG
        if ($type == 1) { # ec number
            my ($ecdesc) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ?",
                                             {}, $objectId);
            unless ($objectId =~ m/-/) {
              $url = "myFitShow.cgi?gene=ec:$objectId";
              $url .= "&orgId=$orgId" if $orgId;
            }
            my $ecobj = { ecnum => $objectId, ecdesc => $ecdesc || $objectId,
                          coord => $coord, url => $url };
            push @ecnums, $ecobj;
        } elsif ($type == 0 || $type == 2) { # compounds or maps
            $url = "keggmap.cgi?mapId=$objectId" if $type == 2;
            push @links, [ $coord, $url ];
        } else {
            ; # ignore type 3 or above, i.e. reactions
        }
    }
    return { mapId => $mapId, mapdesc => $mapdesc,
             links => \@links, ecnums => \@ecnums };
}

sub DrawKEGGMap($) {
    my ($map) = @_;
    my $mapId = $map->{mapId};
    print
        qq{<DIV style="position: relative; left:0; top:0;">},
        qq{<IMG src="../kegg/maps/map$mapId.png" usemap="#imagemap" style="position: relative; top: 0; left: 0;">\n};
    print qq{<MAP name="imagemap" id="imagemap" />\n};
    foreach my $row (@{ $map->{links} }) {
        my ($coord, $url) = @$row;
        my ($coordtype, $pos) = split /:/, $coord;
        next unless $coordtype eq "rect" || $coordtype eq "circ";
        print qq{<AREA SHAPE="$coordtype" coords="$pos" href="$url" />\n};
    }
    print qq{</MAP>};

    foreach my $row (@{ $map->{ecnums}}) {
        my $coord = $row->{coord};
        my ($coordtype, $pos) = split /:/, $coord;
        next unless $coordtype eq "rect";
        my ($left,$top,$right,$bottom) = split /,/, $pos;
        die "Invalid coords $coord" unless defined $top;
        my $width = $right-$left+1;
        my $height = $bottom-$top+1;
        my $alt = $row->{ecdesc};
        my $url = $row->{url};

        my $borderstyle = "";
        $borderstyle = " border: solid red; "
            if $row->{highlight};

        if ($row->{colors}) {
            my @colors = @{ $row->{colors} };
            my $n = scalar(@colors);
            foreach my $i (0..($n-1)) {
                my $x1 = int(0.5 + $left + ($right-$left) * $i/$n);
                my $x2 = int(0.5 + $left + ($right-$left) * ($i+1)/$n);
                my $xw = $x2 - $x1 + 1;
                print qq{<A style="position:absolute; top:${top}px; left:${x1}px; width:${xw}px; height:${height}px; background-color: $colors[$i]; opacity: 0.6; $borderstyle" title="$alt" href="$url"></A>\n};
            }
        } else {
            my $bg = "rgba(255,255,255,0)";
            if ($row->{mask}) {
                $bg = "rgba(0,0,0,0.3)";
            }
            print qq{<A style="position:absolute; top:${top}px; left:${left}px; width:${width}px; height:${height}px; background-color: ${bg}; $borderstyle" title="$alt" href="$url"></A>\n};
        }
    }
    print "</DIV>\n"; # close the div that contains the image
}
