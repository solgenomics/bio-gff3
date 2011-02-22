package Bio::GFF3::Transform::SyncDirectives;
use strict;
use warnings;

use File::Temp;
use File::ReadBackwards;

use Bio::GFF3::LowLevel ();

require Exporter;
our @ISA = ( 'Exporter' );
our @EXPORT_OK = ( 'gff3_add_sync_directives' );

sub gff3_add_sync_directives {
    my ( $out_fh, @files ) = @_;

    my $tempfile = File::Temp->new;

    my %open_parent_rels;
    for my $file ( @files ) {
        my $fh = File::ReadBackwards->new( $file );
        while ( my $line = $fh->readline ) {
            $tempfile->print( $line ) unless $line =~ /^###\s*$/;
            unless( $line =~ /^#/ ) {
                if ( my ( $attr ) = $line =~ / \t ([^\t]+) $/x ) {
                    $attr = Bio::GFF3::LowLevel::gff3_parse_attributes( $attr );
                    if ( $attr->{Parent} ) {
                        for ( @{$attr->{Parent}} ) {
                            $open_parent_rels{ $_ } = 1;
                        }
                    }
                    if ( $attr->{ID} ) {
                        for ( @{$attr->{ID}} ) {
                            delete $open_parent_rels{ $_ };
                        }
                    }
                }
                $tempfile->print( "###\n" ) unless %open_parent_rels;
            }
        }
    }
    $tempfile->close;

    my $temp_backwards = File::ReadBackwards->new( "$tempfile" );
    # print up to and not including the first sync mark (to get rid of the
    # unnecessary first one
    while ( my $line = $temp_backwards->readline ) {
        last if $line =~ /^###$/;
        print $out_fh $line;
    }
    while ( my $line = $temp_backwards->readline ) {
        print $out_fh $line;
    }
}

1;

__END__

=head1 NAME

Bio::GFF3::Transform::SyncDirectives - insert sync (###) directives
into an existing GFF3 file.

=head1 SYNOPSIS

    use Bio::GFF3::Transform::SyncDirectives 'gff3_add_sync_directives';

    my @input_files = ( 'input1.gff3', 'input2.gff3' );
    open my $output_fh, '>', 'myoutputfile.gff3';
    gff3_add_sync_directives( $output_fh, @input_files );

=head1 FUNCTIONS

All functions below are EXPORT_OK.

=head2 gff3_add_sync_directives( $out_filehandle, @files )

Read GFF3 from the given files, add as many sync directives (###) as
possible, and print the resulting GFF3 to the given output filehandle.
Existing sync directives will not be preserved.

=cut

