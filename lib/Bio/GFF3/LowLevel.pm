package Bio::GFF3::LowLevel;
# ABSTRACT: fast, low-level functions for parsing and formatting GFF3

use strict;

use Scalar::Util ();
use URI::Escape ();

=head1 SYNOPSIS

  use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;

  open my $gff3_fh, 'myfile.gff3' or die;
  while( <$gff3_fh> ) {
    next if /^#/;
    my $feat = gff3_parse_feature( $_ );
  }

=head1 DESCRIPTION

These are low-level, fast functions for parsing GFF version 3 files.
All they do is convert back and forth between low-level Perl data
structures and GFF3 text.

Sometimes this is what you need when you are just doing simple
transformations on GFF3.  I found myself writing these functions over
and over again, until I finally got fed up enough to just package them
up properly.

These functions do no validation, do not reconstruct feature
hierarchies, or anything like that.  If you want that, use
L<Bio::FeatureIO>.

All of the functions in this module are EXPORT_OK, meaning that you
can add their name after using this module to make them available in
your namespace.

=cut

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
  gff3_parse_feature
  gff3_parse_attributes
  gff3_format_feature
  gff3_format_attributes
  gff3_escape
  gff3_unescape
);

my @gff3_field_names = qw(
    seq_id
    source
    type
    start
    end
    score
    strand
    phase
    attributes
);

=func gff3_parse_feature( $line )

Given a string containing a GFF3 feature line (i.e. not a comment),
parses it and returns a hashref of its information, of the form:

    {
        seq_id => 'chr02',
        source => 'AUGUSTUS',
        type   => 'transcript',
        start  => '23486',
        end    => '48209',
        score  => '0.02',
        strand => '+',
        phase  => undef,
        attributes => {
            ID => [
                'chr02.g3.t1'
              ],
            Parent => [
                'chr02.g3'
              ],
          },
    }

Note that all values are simple scalars, except for C<attributes>,
which is a hashref as returned by L</gff3_parse_attributes> below.

Unescaping is performed according to the GFF3 specification.

=cut

sub gff3_parse_feature {
  my ( $line ) = @_;
  chomp $line;

  my @f = split /\t/, $line;
  for( 0..8 ) {
      if( defined $f[$_] && $f[$_] eq '.' ) {
          $f[$_] = undef;
      }
  }
  # don't unescape the attr column, that is parsed separately
  for( 0..7 ) {
      $f[$_] = gff3_unescape( $f[$_] );
  }

  $f[8] = gff3_parse_attributes( $f[8] );
  my %parsed;
  @parsed{@gff3_field_names} = @f;
  return \%parsed;
}


=func gff3_parse_attributes( $attr_string )

Given a GFF3 attribute string, parse it and return a hashref of its
data, of the form:

    {
      'attribute_name' => [ value, value, ... ],
      ...
    }

Always returns a hashref.  If the passed attribute string is
undefined, or ".", the hashref returned will be empty.  Attribute
values are always arrayrefs, even if they have only one value.

=cut

sub gff3_parse_attributes {
    my ( $attr_string ) = @_;

    return {} if !defined $attr_string || $attr_string eq '.';

    chomp $attr_string;

    my %attrs;
    for my $a ( split /;/, $attr_string ) {
        next unless $a;
        my ( $name, $values ) = split /=/, $a, 2;
        next unless defined $values;
        push @{$attrs{$name}}, map gff3_unescape($_), split /,/, $values;
    }

    return \%attrs;
}

=func gff3_format_feature( \%fields )

Given a hashref of feature information in the same format returned by
L</gff3_parse_feature> above, constructs a correctly-escaped line of
GFF3 encoding that information.

The line ends with a single newline character, a UNIX-style line
ending, regardless of the local operating system.

=cut

sub gff3_format_feature {
    my ( $f ) = @_;

    my $attr_string = $f->{attributes};
    $attr_string = '.' unless defined $attr_string;

    $attr_string = gff3_format_attributes( $attr_string )
        if    ref( $attr_string ) eq 'HASH'
           && ! Scalar::Util::blessed( $attr_string );

    return join( "\t",
                 ( map { defined $_ ? gff3_escape($_) : '.' }
                   @{$f}{@gff3_field_names[0..7]}
                 ),
                 $attr_string
               )."\n";
}

=func gff3_format_attributes( \%attrs )

Given a hashref of GFF3 attributes in the same format returned by
L</gff3_parse_attributes> above, returns a correctly formatted and
escaped GFF3 attribute string (the 9th column of a GFF3 feature line)
encoding those attributes.

For convenience, single-valued attributes can have simple scalars as
values in the passed hashref.  For example, if a feature has only one
C<ID> attribute (as it should), you can pass C<{ ID =E<gt> 'foo' }>
instead of C<{ ID =E<gt> ['foo'] }}>.

=cut

sub gff3_format_attributes {
  my ( $attr ) = @_;

  return '.' unless defined $attr;

  my $astring = join ';' => (
    map {
      my $key = $_;
      my $val = $attr->{$key};
      $val = [ $val ] unless ref $val;
      "$key=".join( ',', map gff3_escape($_), @$val );
    } sort keys %$attr
 );

  return length $astring ? $astring : '.';
}

=func gff3_escape( $string )

Given a string, escapes special characters in that string according to
the GFF3 specification.

=cut

sub gff3_escape {
    URI::Escape::uri_escape( $_[0], '\n\r\t;=%&,\x00-\x1f\x7f-\xff' )
}

=func gff3_unescape( $string )

Unescapes a GFF3-escaped string.

=cut

sub gff3_unescape {
    URI::Escape::uri_unescape( $_[0] )
}

