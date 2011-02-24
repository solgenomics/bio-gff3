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
    seqid
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

=cut

sub gff3_parse_feature {
  my ( $line ) = @_;
  chomp $line;

  my @f = split /\t/, $line;
  for( 0..8 ) {
      if( $f[$_] eq '.' ) {
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

=cut

sub gff3_parse_attributes {
    my ( $attr_string ) = @_;
    chomp $attr_string;

    return undef if !defined $attr_string || $attr_string eq '.';

    my %attrs;
    for my $a ( split /;/, $attr_string ) {
        next unless $a;
        my ( $name, $values ) = split /=/, $a, 2;
        push @{$attrs{$name}}, map gff3_unescape($_), split /,/, $values;
    }

    return \%attrs;
}

=func gff3_format_feature( \@fields, \%attrs )

=cut

sub gff3_format_feature {
    my ( $f ) = @_;

    my $attr_string = $f->{attributes};
    $attr_string = '.' unless defined $attr_string;

    $attr_string = gff3_format_attributes( $attr_string )
        if ref( $attr_string ) eq 'HASH' && !Scalar::Util::blessed( $attr_string );

    return join( "\t",
                 ( map { defined $_ ? gff3_escape($_) : '.' }
                   @{$f}{@gff3_field_names[0..7]}
                 ),
                 $attr_string
               )."\n";
}

=func gff3_format_attributes( \%attrs )

=cut

sub gff3_format_attributes {
  my ( $attr ) = @_;

  return join ';' => (
    map {
      my $key = $_;
      my $val = $attr->{$key};
      $val = [ $val ] unless ref $val;
      "$key=".join( ',', map gff3_escape($_), @$val );
    } sort keys %$attr
 );

}

=func gff3_escape( $string )

=cut

sub gff3_escape {
    URI::Escape::uri_escape( $_[0], '\n\r\t;=%&,\x00-\x1f\x7f-\xff' )
}

=func gff3_unescape( $string )

=cut

sub gff3_unescape {
    URI::Escape::uri_unescape( $_[0] )
}

