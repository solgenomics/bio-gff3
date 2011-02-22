#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use File::Temp;
use File::ReadBackwards;

use Bio::GFF3::LowLevel 'gff3_parse_attributes';

@ARGV or pod2usage();

my $tempfile = File::Temp->new;

my %open_parent_rels;
for my $file ( @ARGV ) {
    my $fh = File::ReadBackwards->new( $file );
    while( my $line = $fh->readline ) {
        $tempfile->print( $line );
        unless( $line =~ /^#/ ) {
            if( my ( $attr ) = $line =~ / \t ([^\t]+) $/x ) {
                $attr = gff3_parse_attributes( $attr );
                if( $attr->{Parent} ) {
                    $open_parent_rels{ $_ } = 1
                        for @{$attr->{Parent}};
                }
                if( $attr->{ID} ) {
                    delete $open_parent_rels{ $_ }
                        for @{$attr->{ID}};
                }
            }
            $tempfile->print( "###\n" ) unless %open_parent_rels;
        }
    }
}
$tempfile->close;

my $temp_backwards = File::ReadBackwards->new( "$tempfile" );
my $first_sync = 1;
# print up to and not including the first sync mark (to get rid of the
# unnecessary first one
while( my $line = $temp_backwards->readline ) {
    last if $line =~ /^###$/;
    print $line;
}
while( my $line = $temp_backwards->readline ) {
    print $line;
}

__END__

=head1 NAME

gff3_insert_sync_marks.pl - efficiently insert sync (###) marks into a
GFF3 file.  Prints resulting gff3 to STDOUT.

=head1 USAGE

    gff3_insert_sync_marks.pl  file.gff3  anotherfile.gff3 ... > with_syncs.gff3

=cut
