use strict;
use warnings;
use Test::More;

use IO::Scalar;
use File::Temp;

use Bio::GFF3::Transform::SyncDirectives 'gff3_add_sync_directives';

{
    my $test_gff3 = 't/data/gff3_with_syncs.gff3';
    my $t1 = file_without_syncs( $test_gff3 );
    my $out = undef;
    gff3_add_sync_directives( IO::Scalar->new( \$out ), $t1 );

    is( $out, read_file( $test_gff3 ), 'got right sync marks' );
}

done_testing;

sub read_file {
    open my $f, '<', shift or die "$!";
    local $/;
    return <$f>;
}

sub file_without_syncs {
    my $t = File::Temp->new;
    open my $f, '<', +shift or die;
    while( my $line = <$f> ) {
        $t->print( $line ) unless $line =~ /^###$/;
    }
    $t->close;
    return $t;
}
