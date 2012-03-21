package Bio::GFF3::LowLevel::Parser;
# ABSTRACT: a fast, low-level gff3 parser
use strict;
use warnings;
use Carp;

use IO::Handle ();
use Scalar::Util ();

use Bio::GFF3::LowLevel ();

=head1 SYNOPSIS

  my $p = Bio::GFF3::LowLevel::Parser->new( $file_or_fh );
  while( my $i = $p->next_item ) {
      if( exists $i->{seq_id} ) {

      }
      elsif( $i->{directive} ) {

      }
      elsif( $i->{FASTA_fh} ) {

      }
      else {
          die 'this should never happen!';
      }
  }

=head1 DESCRIPTION

This is a fast, low-level parser for Generic Feature Format, version 3
(GFF3).  It is a low-level parser, it only returns dumb hashrefs.  It
B<does> construct feature hierarchies, however, storing an arrayref of
child features under a C<children> key of the parent feature's
hashref.

=cut

=func new

=cut

sub new {
    my $class = shift;
    return bless {

        filethings  => \@_,
        filehandles => [ map $class->_open($_), @_ ],

        item_buffer => [],

        under_construction_top_level_features => [],
        under_construction_by_id              => {},

    }, $class;
}
sub _open {
    my ( $class, $thing ) = @_;
    return $thing if ref $thing eq 'GLOB' || Scalar::Util::blessed( $thing ) && $thing->can('getline');
    open my $f, '<', $thing or croak "$! opening '$thing' for reading";
    return $f;
}

=func next_item

=cut

sub next_item {
    my ( $self ) = @_;
    my $item_buffer = $self->{item_buffer};

    # try to get more items if the buffer is empty
    $self->_buffer_items unless @$item_buffer;

    # return the next item if we have some
    return shift @$item_buffer if @$item_buffer;

    # if we were not able to get any more items, return nothing
    return;
}

## get and parse lines from the files(s) to add at least one item to
## the buffer
sub _buffer_items {
    my ( $self ) = @_;

    my $item_buffer = $self->{item_buffer};

    while( my $line = $self->_next_line ) {
        if( $line =~ /^ \s* [^#\s>] /x ) { #< feature line, most common case
            my $f = Bio::GFF3::LowLevel::gff3_parse_feature( $line );
            $f->{child_features} = [];
            $self->_buffer_feature( $f );
            return if @$item_buffer; #< return if we were able to buffer the feature for returning
        }
        # directive or comment
        elsif( my ( $hashsigns, $contents ) = $line =~ /^ \s* (\#+) (.*) /x ) {
            if( length $hashsigns == 3 ) { #< sync directive, all forward-references are resolved.
                $self->_buffer_all_under_construction_features;
            }
            elsif( length $hashsigns == 2 ) {
                my $directive = Bio::GFF3::LowLevel::gff3_parse_directive( $line );
                if( $directive->{directive} eq 'FASTA' ) {
                    $self->_buffer_all_under_construction_features;
                    push @$item_buffer, { directive => 'FASTA', filehandle => shift @{$self->{filehandles} } };
                    shift @{$self->{filethings}};
                } else {
                    push @$item_buffer, $directive;
                }
            }
            else {
                $contents =~ s/\s*$//;
                push @$item_buffer, { comment => $contents };
            }
            return;
        }
        elsif( $line =~ /^ \s* $/x ) {
            # blank line, do nothing
        }
        elsif( $line =~ /^ \s* > /x ) {
            # implicit beginning of a FASTA section.  a very stupid
            # idea to include this in the format spec.  increases
            # implementation complexity by a lot.
            $self->_buffer_all_under_construction_features;
            push @$item_buffer, $self->_handle_implicit_fasta_start( $line );
        }
        else { # it's a parse error
            chomp $line;
            croak "$self->{filethings}[0]:$.: parse error.  Cannot parse '$line'.";
        }
    }

    # if we are out of lines, buffer all under-construction features
    $self->_buffer_all_under_construction_features;
}

## take all under-construction features and put them in the
## item_buffer to be output
sub _buffer_all_under_construction_features {
    my ( $self ) = @_;

    push @{$self->{item_buffer}}, @{$self->{under_construction_top_level_features}};

    $self->{under_construction_top_level_features} = [];
    $self->{under_construction_by_id}     = {};
}


## get the next line from our file(s), returning nothing if we are out
## of lines and files
sub _next_line {
    my ( $self ) = @_;
    my $filehandles = $self->{filehandles};
    while( @$filehandles ) {
        my $line = $filehandles->[0]->getline;
        return $line if $line;
        shift @$filehandles;
        shift @{$self->{filethings}};
    }
    return;
}

## do the right thing with a newly-parsed feature line
sub _buffer_feature {
    my ( $self, $feature ) = @_;

    my $ids     = $feature->{attributes}{ID}     || [];
    my $parents = $feature->{attributes}{Parent} || [];

    if( !@$ids && !@$parents ) {
        # if it has no IDs or Parents, it's independent, so just
        # put it in the output buffer
        push @{$self->{item_buffer}}, $feature;
        return;
    }

    if( @$ids ) {
        if( !@$parents ) {
            push @{ $self->{under_construction_top_level_features} }, $feature;
        }
        for my $id ( @$ids ) {
            $self->{under_construction_by_id}{$id} = $feature;
        }
    }

    for my $parent_id ( @$parents ) {
        my $parent = $self->{under_construction_by_id}{$parent_id}
            or die( "No feature found with ID=$parent_id" );

        push @{$parent->{child_features}}, $feature;
    }
}

sub _handle_implicit_fasta_start {
    my ( $self, $line ) = @_;
    require POSIX;
    require IO::Pipe;
    my $pipe = IO::Pipe->new;
    unless( fork ) {
        $pipe->writer;
        my $fh = $self->{filehandles}[0];
        undef $self;
        $pipe->print($line);
        while( $line = $fh->getline ) {
            $pipe->print( $line );
        }
        $pipe->close;
        POSIX::_exit(0);
    }
    $pipe->reader;
    shift @$_ for $self->{filehandles}, $self->{filethings};
    return { directive => 'FASTA', filehandle => $pipe };
}

1;
