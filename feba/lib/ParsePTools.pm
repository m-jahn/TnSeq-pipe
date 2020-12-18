# Parse PathwayTools flat files
package ParsePTools;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(ParsePTools);

# Read a record from a pathway tools attribute-value file
# and returns the result as a reference to a hash of
# attribute_name => list of { value => value, annotation_name => anno }
# where the key 'value' is the main attribute value and the annotation names vary.
#
# Multiple values for an attribute are supported, but if there are
# multiple annotations with the same name for an attribute/value pair,
# then only the first one is reported.
#
# Parsing is based on documentation from
# https://bioinformatics.ai.sri.com/ptools/flatfile-format.html
sub ParsePTools {
  my ($fh) = shift;
  die unless defined $fh;
  local($/) = "\n//\n";
  my $record = <$fh>;
  return undef if !defined $record;
  my @lines = split /\n/, $record;
  foreach (@lines) { chomp; }
  # Remove comments
  @lines = grep { ! m/^#/ } @lines;
  return undef if @lines == 0;
  $lines[-1] eq "//" || die "Last lines in record does not match //";
  pop @lines;
  # Merge trailing lines
  my @merged = ();
  foreach my $line (@lines) {
    if ($line =~ m!^/!) {
      die "Continuation line with no preceding line: $line"
        unless @merged > 0;
      $merged[-1] .= " " . substr($line, 1);
    } else {
      push @merged, $line;
    }
  }
  my $last_attr = undef;
  my $out = {};
  foreach my $line (@merged) {
    $line =~ m/^(\S+) - (.*)$/ || die "Cannot parse attribute or annotation from $line";
    my ($attr, $value) = ($1,$2);
    if ($attr =~ m/^\^/) {
      my $anno = substr($attr, 1);
      die "Annotation with no preceding attribute: $line" unless defined $last_attr;
      my $h = $out->{$last_attr}[-1];
      if (exists $h->{$anno}) {
        # print STDERR "Duplicate annotation for $last_attr $anno\n";
        ;
      } else {
        $h->{$anno} = $value;
      }
    } else {
      push @{ $out->{$attr} }, { 'value' => $value };
      $last_attr = $attr;
    }
  }
  return $out;
}
