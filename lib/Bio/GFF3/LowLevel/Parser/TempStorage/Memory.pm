package Bio::GFF3::LowLevel::Parser::TempStorage::Memory;
use strict;
use warnings;

sub new {
    bless {
        # features that are ready to go out and be flushed
        item_buffer => [],

        # features that we have to keep on hand for now because they
        # might be referenced by something else
        under_construction_top_level => [],
        # index of the above by ID
        under_construction_by_id => {},

        # features that reference something we have not seen yet
        # structured as:
        # {  'some_id' => {
        #     'Parent' => [ orphans that have a Parent attr referencing it ],
        #     'Derives_from' => [ orphans that have a Derives_from attr referencing it ],
        # }
        under_construction_orphans => {},

    }, shift
}

sub output_buffer_add {
    push @{$_[0]->{item_buffer}}, $_[1];
}
sub output_buffer_size {
    scalar @{ $_[0]->{item_buffer} }
}
sub output_buffer_next {
    shift @{ $_[0]->{item_buffer} };
}


sub flush_under_construction {
    my ( $self ) = @_;

    # since the under_construction_top_level buffer is likely to be
    # much larger than the item_buffer, we swap them and unshift the
    # existing buffer onto it to avoid a big copy.
    my $old_buffer = $self->{item_buffer};
    $self->{item_buffer} = $self->{under_construction_top_level};
    unshift @{$self->{item_buffer}}, @$old_buffer;
    undef $old_buffer;

    $self->{under_construction_top_level} = [];
    $self->{under_construction_by_id} = {};

    # if we have any orphans hanging around still, this is a problem. die with a parse error
    if( grep %$_, values %{$self->{under_construction_orphans}} ) {
        require Data::Dumper; local $Data::Dumper::Terse = 1;
        die "parse error: orphans ", Data::Dumper::Dumper($self->{under_construction_orphans}); # TODO: make this better
    }
}
sub under_construction_get {
    $_[0]->{under_construction_by_id}{$_[1]}
}

sub under_construction_add {
    my ( $self, $feature, $id, $is_top_level ) = @_;
    if( $is_top_level ) {
        push @{ $self->{under_construction_top_level} }, $feature;
    }
    $self->{under_construction_by_id}{$id} = $feature;
}

sub orphans_get {
    return $_[0]->{under_construction_orphans}{$_[1]};
}

sub orphans_add {
    my ( $self, $to_id, $attrname, $feature ) = @_;
    push @{ $self->{under_construction_orphans}{$to_id}{$attrname} ||= [] }, $feature;
}

1;
