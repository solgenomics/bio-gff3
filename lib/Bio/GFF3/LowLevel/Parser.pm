package Bio::GFF3::LowLevel::Parser;
# ABSTRACT: a fast, low-level gff3 parser

use strict;
use warnings;
use Carp;

use IO::Handle ();
use Scalar::Util ();

use List::MoreUtils ();

use Bio::GFF3::LowLevel ();
use Bio::GFF3::LowLevel::Parser::TempStorage::Memory ();
use Bio::GFF3::LowLevel::Parser::TempStorage::DBM ();

=head1 SYNOPSIS

  my $p = Bio::GFF3::LowLevel::Parser->open( $file_or_fh );

  while( my $i = $p->next_item ) {

      if( ref $i eq 'ARRAY' ) {
          ## $i is an arrayref of feature lines that have the same ID,
          ## in the same format as returned by
          ## Bio::GFF3::LowLevel::gff3_parse_feature
          for my $f (@$i) {
             # for each location of this feature
             # do something with it
          }
      }
      elsif( $i->{directive} ) {
          if( $i->{directive} eq 'FASTA' ) {
              my $fasta_filehandle = $i->{filehandle};
              ## parse the FASTA in the filehandle with BioPerl or
              ## however you want.  or ignore it.
          }
          elsif( $i->{directive} eq 'gff-version' ) {
              print "it says it is GFF version $i->{value}\n";
          }
          elsif( $i->{directive} eq 'sequence-region' ) {
              print( "found a sequence-region, sequence $i->{seq_id},",
                     " from $i->{start} to $i->{end}\n"
                   );
          }
      }
      elsif( $i->{comment} ) {
          ## this is a comment in your GFF3 file, in case you want to do
          ## something with it.
          print "that comment said: '$i->{comment}'\n";
      }
      else {
          die 'this should never happen!';
      }

  }

=head1 DESCRIPTION

This is a fast, low-level parser for Generic Feature Format, version 3
(GFF3).  It is a low-level parser, it only returns dumb hashrefs.  It
B<does> reconstruct feature hierarchies, however, using features'
C<ID>, C<Parent>, and C<Derives_from> attributes, and it B<does> group
together lines with the same ID (i.e. features that have multiple
locations).

=head3 Features

Features are returned as arrayrefs containing one or more (never zero)
feature lines parsed in the same format as
L<Bio::GFF3::LowLevel/gff3_parse_feature>.  Each has some additional
keys for related features: C<child_features> and C<derived_features>,
each of which is a (possibly empty) arrayref of features
(i.e. arrayrefs) that refer to this one as a C<Parent> or claim that
they C<Derives_from> it.

Note that, to make code that uses this parser easier to write, B<all>
features have C<child_features> and C<derived_features> arrayrefs.
This means you don't have to check for the existence of these before
seeing if they have anything in them.

=head3 Directives

Directives are returned as hashrefs, in the same format as
L<Bio::GFF3::LowLevel/gff3_parse_directive>.

=head3 Comments

Comments are parsed into a hashref of the form:

  { comment => 'text of the comment, not including the hash mark(s) and ending newline' }

=cut

=func open( $file_or_filehandle, ... )

Make a new parser object that will parse the GFF3 from all of the files
or filehandles that you give it, as if they were all a single stream.

=cut

sub open {
    my $class = shift;
    return bless {

        filethings  => \@_,
        filehandles => [ map $class->_open($_), @_ ],

        temp => Bio::GFF3::LowLevel::Parser::TempStorage::DBM->new,

        completed_references => {},
    }, $class;
}
sub _open {
    my ( $class, $thing ) = @_;
    return $thing if ref $thing eq 'GLOB' || Scalar::Util::blessed( $thing ) && $thing->can('getline');
    CORE::open my $f, '<', $thing or croak "$! opening '$thing' for reading";
    return $f;
}

=func new

Returns a wrapped copy of this parser that returns data that is backward-compatible with what the 1.0 version of this parser returned.  Do not use in new code.

=cut

sub new {
    my $class = shift;
    require Bio::GFF3::LowLevel::Parser::1_0_backcompat;
    return Bio::GFF3::LowLevel::Parser::1_0_backcompat->new( @_ );
}

=func next_item()

Iterate through all of the items (features, directives, and comments)
in the file(s) given to the parser.  Each item is a returned as a
hashref.

=cut

sub next_item {
    my ( $self ) = @_;

    # try to get more items if the buffer is empty
    $self->_buffer_items unless $self->_buffered_items_count;

    # return the next item if we have some
    return $self->_shift_item if $self->_buffered_items_count;

    # if we were not able to get any more items, return nothing
    return;
}

sub _shift_item {
    $_[0]->{temp}->output_buffer_next;
}
sub _buffer_item {
    $_[0]->{temp}->output_buffer_add( $_[1] );
}
sub _buffered_items_count {
    $_[0]->{temp}->output_buffer_size;
}


## get and parse lines from the files(s) to add at least one item to
## the buffer
sub _buffer_items {
    my ( $self ) = @_;

    while( my $line = $self->_next_line ) {
        if( $line =~ /^ \s* [^#\s>] /x ) { #< feature line, most common case
            my $f = Bio::GFF3::LowLevel::gff3_parse_feature( $line );
            $self->_buffer_feature( $f );
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
                    $self->_buffer_item( { directive => 'FASTA', filehandle => shift @{$self->{filehandles} } });
                    shift @{$self->{filethings}};
                } else {
                    $self->_buffer_item( $directive );
                }
            }
            else {
                $contents =~ s/\s*$//;
                $self->_buffer_item( { comment => $contents } );
            }
        }
        elsif( $line =~ /^ \s* $/x ) {
            # blank line, do nothing
        }
        elsif( $line =~ /^ \s* > /x ) {
            # implicit beginning of a FASTA section.  a very stupid
            # idea to include this in the format spec.  increases
            # implementation complexity by a lot.
            $self->_buffer_all_under_construction_features;
            $self->_buffer_item( $self->_handle_implicit_fasta_start( $line ) );
        }
        else { # it's a parse error
            chomp $line;
            croak "$self->{filethings}[0]:$.: parse error.  Cannot parse '$line'.";
        }

        # return now if we were able to find some things to put in the
        # output buffer
        return if $self->_buffered_items_count;
    }

    # if we are out of lines, buffer all under-construction features
    $self->_buffer_all_under_construction_features;
}

## take all under-construction features and put them in the
## item_buffer to be output
sub _buffer_all_under_construction_features {
    my ( $self ) = @_;
    $self->{temp}->flush_under_construction;
    $self->{completed_references} = {};
}

## get the next line from our file(s), returning nothing if we are out
## of lines and files
sub _next_line {
    no warnings;
    # fast code path for reading a line from the first filehandle,
    my $first_fh = $_[0]->{filehandles}[0];
    return <$first_fh> || do {
        # slower case where we are at the end, or need to change
        # filehandles
        my ( $self ) = @_;
        my $filehandles = $self->{filehandles};
        while ( @$filehandles ) {
            my $line = $filehandles->[0]->getline;
            return $line if $line;
            shift @$filehandles;
            shift @{$self->{filethings}};
        }
        return;
    };
}

my %container_attributes = ( Parent => 'child_features', Derives_from => 'derived_features' );
## do the right thing with a newly-parsed feature line
sub _buffer_feature {
    my ( $self, $feature_line ) = @_;

    $feature_line->{'child_features'}   = [];
    $feature_line->{'derived_features'} = [];

    # NOTE: a feature is an arrayref of one or more feature lines.
    my $ids     = $feature_line->{attributes}{ID}     || [];
    my $parents = $feature_line->{attributes}{Parent} || [];
    my $derives = $feature_line->{attributes}{Derives_from} || [];

    if( !@$ids && !@$parents && !@$derives ) {
        # if it has no IDs and does not refer to anything, we can just
        # output it
        $self->_buffer_item( $feature_line );
        return;
    }

    my $feature;
    for my $id ( @$ids ) {
        if( my $existing = $self->{temp}->under_construction_get($id) ) {
            # another location of the same feature
            push @$existing, $feature_line;
            $feature = $existing;
            $self->{temp}->under_construction_update( $id, $feature );
        }
        else {
            # haven't seen it yet
            $feature = [ $feature_line ];
            $self->{temp}->under_construction_add( $feature, $id, !@$parents && !@$derives );

            # see if we have anything buffered that refers to it
            $self->_resolve_references_to( $feature, $id );
        }
    }

    # try to resolve all its references
    $self->_resolve_references_from( $feature || [ $feature_line ], { Parent => $parents, Derives_from => $derives }, $ids );
}

sub _resolve_references_to {
    my ( $self, $feature, $id ) = @_;
    my $references = $self->{temp}->orphans_get( $id )
        or return;
    for my $attrname ( keys %$references ) {
        my $pname = $container_attributes{$attrname} || lc $attrname;
        for my $loc ( @$feature ) {
            push @{ $loc->{$pname} },
                  @{ delete $references->{$attrname} };
        }
    }
    $self->{temp}->under_construction_update( $id, $feature );
}
sub _resolve_references_from {
    my ( $self, $feature, $references, $ids ) = @_;
    # go through our references
    #  if we have the feature under construction, put this feature in the right place
    #  otherwise, put this feature in the right slot in the orphans

    for my $attrname ( keys %$references ) {
        my $pname;
        for my $to_id ( @{ $references->{ $attrname } } ) {
            if( my $other_feature = $self->{temp}->under_construction_get( $to_id ) ) {
                $pname ||= $container_attributes{$attrname} || lc $attrname;
                unless( grep $self->{completed_references}{$_}{$attrname}{$to_id}++, @$ids ) {
                    for my $loc ( @$other_feature ) {
                        push @{ $loc->{ $pname } }, $feature;
                    }
                }
                $self->{temp}->under_construction_update( $to_id, $other_feature );
            }
            else {
                $self->{temp}->orphans_add( $to_id, $attrname, $feature );
            }
        }
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
