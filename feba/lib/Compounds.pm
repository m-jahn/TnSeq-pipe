# Compounds.pm -- utilities for reading the media and compounds metadata from (by default)
# feba/metadata/Compounds.tsv, feba/metadata/media, feba/metadata/mixes
#
# It maintains a lot of state in global variables.

package Compounds;
require Exporter;
use strict;
use FEBA_Utils; # for ReadTable(), ReadColumnNames()

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(LoadCompounds
             GetCompoundList
             FindCompound GetCompoundCAS GetCompoundMW
             ListValidUnits
             LoadMedia
             GetMedias GetMediaComponents GetUndefMedia GetUnknownComponents WarnReusedComponents
             GetMixes GetMixComponents
             SynToKey GetSynonymMap GetCompoundFormatDoc);

sub LoadCompounds($); # metadata directory; no return value
sub FindCompound($); # compound or synonym => compound or undef
sub GetCompoundCAS($); # compound => CAS or ""
sub GetCompoundMW($); # compound => MW or ""
sub GetCompoundList(); # returns the list of compounds

sub LoadMedia($); # metadata directory; no return value
sub GetMedias(); # returns list of media names that have definitions
sub GetUndefMedia(); # returns list of media names that lack definitions
sub GetMediaAttributes($); # media name => hash of attributes such as Description or Minimal
sub GetMediaComponents($); # for each component, a list of [ compound, concentration, units, mix ]
#	where mix indicates which mix it was from (if any)
sub WarnReusedComponents(); # report to STDOUT on components that are in a medium more than once

sub GetMixComponents($); # for a given mix, a list of [ compounds, concentration, units ]


# From a compound name or synonym to a key to look up.
# Only alphanumeric characters, +, or - are considered -- everything else is removed. Also, case is ignored.
# Ignoring whitespace and case means that many variant names for compounds can be handled without introducing more synonyms
sub SynToKey($);

# Local variables for compound information
my %compounds = (); # compound to list of [ compound name, cas, MW ]; cas or MW are "" (not undef) if unknown
my ($I_COMPOUND,$I_CAS,$I_MW) = 0..2;
my %synonyms = (); # processed synonym (from SynToKey) => compound name
sub GetSynonymMap() { return \%synonyms; }

# Local variables for media information
my %media = (); # media => list of [ compound_id, number, units ]
my %mediaAttr = (); # media => Description or Minimal => value
my %mediaExclude = (); # media => excluded compound => 1
my %mix = (); # mix => list of [ compound_id, number, units ]
my %mixAttr = (); # mix => Description or X => value
my %unknownComponents = (); # media components with no match in the compounds table
my @undefMedia = (); # media with no definition
my %reuseComponents = (); # component => media => 1 if it is reused

sub LoadCompounds($) {
    my ($compoundsDir) = @_;
    my $compoundsFile = "$compoundsDir/Compounds.tsv";
    my @req = qw{Compound CAS FW Synonyms};
    my @compounds = &ReadTable($compoundsFile, @req);
    my @headers = &ReadColumnNames($compoundsFile);
    die "No rows in $compoundsFile" unless @compounds > 0;

    foreach my $row (@compounds) {
        my $compound = $row->{"Compound"};
        die "Duplicate compound id $compound" if exists $compounds{$compound};
        my $mw = $row->{"FW"};
        $mw = "" if $mw eq "NA";
        die "Invalid weight $mw for compound $compound in $compoundsFile\n"
            unless $mw eq "" || $mw =~ m/^[0-9]+[.]?[0-9]*$/;
        my $cas = $row->{"CAS"};
        $cas =~ s/ +$//;
        $cas =~ s/^ +//;
        $cas = "" if $cas eq "NA";
        $cas =~ s! / .*!!; # occasionally have more than one CAS number separated by a /, just keep the first one
        die "Invalid cas number '$cas' for compound $compound in $compoundsFile\n"
            unless $cas eq "" || $cas =~ m/^\d+[0-9-]*\d$/;
        $compounds{$compound} = [ $compound, $cas, $mw ];

        my @syns = $compound;
        if ($row->{Synonyms} =~ m/;/) {
          # allow semicolon separators instead of comma separators
          push @syns, split /; /, $row->{Synonyms};
        } else {
          push @syns, split /, /, $row->{Synonyms};
        }
        foreach my $syn (@syns) {
            my $key = SynToKey($syn);
            next if $key eq "";
            print "Warning: non-unique synonym $syn: $compound or $synonyms{$key}\n"
                if exists $synonyms{$key} && $synonyms{$key} ne $compound;
            $synonyms{$key} = $compound;
        }
    }
}

sub FindCompound($) {
    my ($syn) = @_;
    return $syn if exists $compounds{$syn};
    my $key = SynToKey($syn);
    return $synonyms{$key} if exists $synonyms{$key};
    return undef;
}

sub GetCompoundCAS($) {
    my ($compound) = @_;
    return exists $compounds{$compound} ? $compounds{$compound}[$I_CAS] : "";
}

sub GetCompoundMW($) {
    my ($compound) = @_;
    return exists $compounds{$compound} ? $compounds{$compound}[$I_MW] : "";
}

sub GetCompoundList() {
    return sort keys %compounds;
}

sub SynToKey($) {
    my ($syn) = @_;
    $syn =~ s/[^a-zA-Z0-9+-]//g;
    return(lc($syn));
}

### Media functions

# Returns two hashes:
# media => list of components
# media => attributes
# Handles both the media file and the mixes file
# Does NOT look up compound synonyms or otherwise check the results
sub LoadMedia($) {
    my ($metadir) = @_;
    my $mediaFile = "$metadir/media";
    my $mixFile = "$metadir/mixes";
    die "No such file: $mediaFile\n" unless -e $mediaFile;
    die "No such file: $mixFile\n" unless -e $mixFile;
    my ($mediaC,$mediaA) = ParseMediaFile($mediaFile);
    my ($mixC,$mixA) = ParseMediaFile($mixFile);
    %media = %$mediaC;
    %mediaAttr = %$mediaA;
    %mix = %$mixC;
    %mixAttr = %$mixA;

    # Validation:
    # Each media must have a Description
    # media should not also be mixes
    while (my ($media, $attr) = each %mediaAttr) {
        die "No Description for media $media" unless exists $attr->{Description} && $attr->{Description} ne "";
        die "Media $media is also a mix" if exists $mix{$media};
    }
    # Each mix must have a Description and a numeric X value
    # mixes should not also be media
    while (my ($mix, $attr) = each %mixAttr) {
        die "No Description for mix $mix" unless exists $attr->{Description} && $attr->{Description} ne "";
        die "Invalid X for mix $mix" unless exists $attr->{X} && $attr->{X} =~ m/^[0-9]+[.]?[0-9]*$/;
        die "Mix $mix is also a media" if exists $media{$mix};
    }

    # Expand media in terms of other media
    my $nMaxCycles = 100; # not sure if complex cycles will always be detected via substitution
    for (my $nCycle = 0; ; $nCycle++) {
      die "Too many cycles of media expansion -- is there a cycle?\n" if $nCycle >= $nMaxCycles;
      my $nChanges = 0;
      foreach my $media (keys %media) {
        $nChanges += ExpandMedia($media);
      }
      last if $nChanges == 0;
    }

    # Replace compound synonyms with compounds, and record any that are not known or are duplicates
    while (my ($media, $list) = each %media) {
        SetupComponentList($media, $list);
    }
    while (my ($mix, $list) = each %mix) {
        SetupComponentList($mix, $list);
    }
    # Replace compound synonyms with compounds
    # Warn if excluded compound is in the media (and remove it)
    # Define "minus" mixes as needed
    while (my ($media, $hash) = each %mediaExclude) {
      SetupExclude($media, $hash);
    }
}

sub SetupComponentList($$) {
    my ($media, $list) = @_;
    my $isMedia = exists $media{$media};
    my $isMix = exists $mix{$media};
    die "Unknown $media" unless $isMedia || $isMix;

    my $COMPOUND = 0;
    # transfer synonyms or record that it is unknown
    foreach my $row (@$list) {
        my ($orig,$undef,$units) = @$row;
        if ($isMedia && exists $mix{$orig}) {
            # leave as is, but check that X is specified
            die "Mix must be included with X units" unless $units eq "X";
        } else {
            my $compound = FindCompound($orig);
            if (defined $compound) {
                $row->[$COMPOUND] = $compound;
            } else {
                $unknownComponents{$row->[$COMPOUND]} = 1;
            }
        }
    }
    # record repeat entries
    my %seen = ();
    foreach my $row (@$list) {
        my $compound = $row->[$COMPOUND];
        $reuseComponents{$compound}{$media} = 1 if exists $seen{$compound};
        $seen{$compound} = 1;
    }
}

sub SetupExclude($$) {
  my ($media, $excludeHash) = @_;
  my @excluded = ();
  foreach my $orig (keys %$excludeHash) {
    my $compound = FindCompound($orig);
    if (defined $compound) {
      push @excluded, $compound;
    } else {
      $unknownComponents{$orig} = 1;
    }
  }
  # hash of excluded compound => number of alterations
  my %excluded = map { $_ => 0 } @excluded;

  # And update the media and any incorporated mixes to actually exclude
  my @updatedComponents = ();
  foreach my $component (@{ $media{$media} }) {
    my ($compound, $number, $units) = @$component;
    if (exists $excluded{$compound}) {
      $excluded{$compound}++; # exclude it
    } elsif (exists $mix{$compound}) {
      # $compound is a mix -- check if it needs to be altered
      my %mixExcluded = ();
      my @mixKeep = ();
      foreach my $mixComponent (@{ $mix{$compound} }) {
        my ($mixCompound, undef, undef) = @$mixComponent;
        if (exists $excluded{$mixCompound}) {
          $mixExcluded{$mixCompound} = 1;
        } else {
          push @mixKeep, $mixComponent;
        }
      }
      if (scalar(keys %mixExcluded) > 0) {
        # Define the new mix and use it instead
        my $minusString = join(" ", map { "minus $_" } sort keys %mixExcluded);
        my $mixNew = "$compound $minusString";
        if (exists $mix{$mixNew}) {
          die "Wrong number of components for mix $mixNew which can also be built via exclude from $compound"
            unless scalar(@{ $mix{$mixNew} }) == scalar(@mixKeep);
        } else {
          $mix{$mixNew} = \@mixKeep;
          $mixAttr{$mixNew} = $mixAttr{$compound};
          $mixAttr{$mixNew}{Description} .= " $minusString";
        }
        push @updatedComponents, [ $mixNew, $number, $units ];
        foreach my $mixCompound (keys %mixExcluded) {
          $excluded{$mixCompound}++;
        }
      } else {
        # Keep the mix as is
        push @updatedComponents, $component;
      }
    } else {
      # not excluded or a mix, keep it
      push @updatedComponents, $component;
    }
  }
  $media{$media} = \@updatedComponents;
  while (my ($compound, $count) = each %excluded) {
    die "$media $compound" unless defined $count;
    print STDERR "Warning: excluding $compound from media $media had no effect\n"
      unless $count > 0;
  }
}

sub ParseMediaFile($) {
    my ($mediaFile) = @_;
    my %comp = ();
    my %attr = ();
    open(MEDIA, "<", $mediaFile) || die "Cannot read $mediaFile";
    my ($COMPOUND,$NUMBER,$UNITS) = 0..2;

    my %validUnits = map { $_ => 1 } ListValidUnits();
    my %validAttr = map { $_ => 1 } qw{Description Minimal X};

    my $curMedia = undef;
    my $readingCompounds = 0;
    while(my $line = <MEDIA>) {
        $line =~ s/[\r\n]+$//; # handle DOS mode files
        $line =~ s/\t+$//; # strip trailing fields that are empty (note this means units *must* be present)
        my @F = split /\t/, $line;

        if (scalar(@F) == 0) {
            $curMedia = undef; # blank lines end media descriptions
        } elsif ($F[0] =~ m/^#/) {
            # skip comment line
            ;
        } elsif (scalar(@F) == 2) {
            my ($attr,$value) = @F;
            if ($attr eq "Media") {
                $curMedia = $F[1];
                $curMedia =~ s/ +$//;
                die "Duplicate media entry for $curMedia" if exists $comp{$curMedia};
                $comp{$curMedia} = [];
                $attr{$curMedia} = {};
                $readingCompounds = 0;
            } elsif (exists $validAttr{$attr}) {
                die "No media id yet at line:\n$line\nin $mediaFile" if !defined $curMedia;
                die "Duplicate attr $attr for media $curMedia" if exists $attr{$curMedia}{$attr};
                $attr{$curMedia}{$attr} = $value;
            } else {
                die "Invalid media attribute $F[0]";
            }
        } elsif (scalar(@F) == 3) {
            die "No media id yet at line:\n$line\nin $mediaFile" if !defined $curMedia;
            if ($F[0] =~ m/^Controlled/ && $F[1] eq "Concentration" && $F[2] eq "Units") {
                $readingCompounds = 1;
            } else {
                die "No compounds header for $curMedia at\n$line\n..." unless $readingCompounds == 1;
                my ($compound, $concentration, $units) = @F;
                $compound =~ s/ +$//; # remove trailing spaces
                $concentration =~ s/^ +//;
                $concentration =~ s/ +$//;
                if ($concentration eq "-" && $units eq "-") {
                  $mediaExclude{$curMedia}{ $compound } = 1;
                } else {
                  $concentration eq "" || $concentration =~ m/^\d+$/
                    || $concentration =~ m/^\d+[.]\d*$/
                    || $concentration =~ m/^\d+[.]?\d*[eE][+-]\d+$/
                    || die "Invalid concentration $concentration in line\n$line\nfor $curMedia\nin $mediaFile";
                  $units =~ s/^ +//;
                  $units =~ s/ +$//;
                  die $line if $units eq "";
                  die "Invalid unit $units for\n$line\nin $curMedia, $mediaFile"
                    unless exists $validUnits{$units};
                  push @{ $comp{$curMedia} }, [ $compound, $concentration, $units ];
                }
            }
        } else {
            die "Wrong number of fields in\n$line\n";
        }
    }
    close(MEDIA) || die "Error reading $mediaFile";
    return (\%comp, \%attr);
}

sub GetMedias() {
    return sort keys %media;
}

# returns undef for unknown media; otherwise, a reference to a list of [compound,concentration,units,mix]
# where mix is empty (not undefined) unless the compound was included indirectly via a mix
sub GetMediaComponents($) {
    my ($media) = @_;
    return undef if !exists $media{$media};
    my $out = [];
    foreach my $row (@{ $media{$media} }) {
        my ($comp,$conc,$units) = @$row;
        if (exists $mix{$comp}) {
            die "Units for mix $comp in media $media are not X" unless $units eq "X";
            die "No X value for mix $comp" unless exists $mixAttr{$comp}{"X"};
            my $rel = $conc / $mixAttr{$comp}{"X"};
            die "Invalid relative X $rel for $comp in $media" unless $rel > 0 && $rel < 1e4;
            foreach my $row2 (@{ $mix{$comp} }) {
                my ($comp2, $conc2, $units2) = @$row2;
                push @$out, [ $comp2, $conc2 * $rel, $units2, $comp ];
            }
        } else {
            push @$out, [ $comp, $conc, $units, "" ];
        }
    }
    return $out;
}

# return number of changes made
sub ExpandMedia($$) {
  my ($media) = @_;
  my @components = ();
  my $nExpand = 0;
  foreach my $row (@{ $media{$media} }) {
    my ($comp,$conc,$units) = @$row;
    die "Invalid definition of medium $media -- it includes itself\n" if $comp eq $media;
    if (exists $media{$comp}) { # $comp is another media
      die "Invalid definition of medium $media -- it includes medium $comp with an excluded compound\n"
        if exists $mediaExclude{$comp};
      die "Invalid definition of medium $media -- it includes medium $comp but units are not X\n"
        unless $units eq "X";
      foreach my $subComponent (@{ $media{$comp} }) {
        my ($comp2,$conc2,$units2) = @$subComponent;
        $conc2 *= $conc; # multiply by X
        push @components, [ $comp2, $conc2, $units2 ];
      }
      $nExpand++;
    } else {
      push @components, $row;
    }
  }
  $media{$media} =  \@components;
  return $nExpand;
}

sub GetMixes() {
    return sort keys %mix;
}

# Like GetMediaComponents but for mixes; returns a list of [compound,concentration,units]
sub GetMixComponents($) {
  my ($mix) = @_;
  return undef if !exists $mix{$mix};
  return $mix{$mix};
}

sub GetUnknownComponents() {
    return sort keys %unknownComponents;
}

sub GetUndefMedia() {
    return @undefMedia;
}

sub WarnReusedComponents() {
    foreach my $compound (sort keys %reuseComponents) {
        foreach my $media (keys %{ $reuseComponents{$compound} }) {
            my @components = grep { $_->[0] eq $compound } @{ $media{$media} };
            die "Not dup compound $compound media $media: " . scalar(@components) if scalar(@components) < 2;
            my @show = map { $_->[1] . " " . $_->[2] } @components;
            print "Multiple use of $compound in $media: " . join(", ", @show)."\n";
        }
    }
}

sub ListValidUnits() {
    return qw{g/L mg/L ug/L M mM uM vol% ml/L X};
}

sub GetCompoundFormatDoc() {
  my $unitString = join(", ", ListValidUnits);
  return <<END
Compounds.tsv file should be tab-delimited with the first five fields
being a unique id, CAS no, source, catalog no, molecular
weight. It should also  contain a field named Synonyms.

Each medium's definition begins like this:
--
Media	media identifier
Description	description for this media
Minimal	TRUE
Controlled vocabulary	Concentration	Units
--
and will then have lines of the form
compound	concentration	units
where valid units are $unitString.

The mixes file is in a similar format  but should begin like this:
--
Media	mix identifier
Description	description for this mix
X	100
Controlled vocabulary	Concentration	Units
--
where X=100 means that the concentration is 100X higher than it would
usually appear in the final media (where mixes are usually used at 1X).

Media can incorporate mixes or other media with the concentration
given as X. (Media are always assumed to be defined as X=1.) Mixes may
not incorporate other media or mixes.

In the list of media components, concentration="-" and units="-" means
to exclude that compound from the media, including from any media or
mixes that the media is built out of. When listing the components of a
media whose mix was modified by exclusion, the mix's name is modified
to "mix minus compound1 minus compound2".

A media with the exclusion feature may not be incorporated into other
media.
END
    ;
}

1;
