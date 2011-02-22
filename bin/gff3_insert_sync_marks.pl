#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Pod::Usage;

use Bio::GFF3::Transform::SyncDirectives 'gff3_add_sync_directives';

@ARGV or pod2usage();

gff3_add_sync_directives( \*STDOUT, @ARGV );

__END__

=head1 NAME

gff3_insert_sync_marks.pl - efficiently insert sync (###) marks into a
GFF3 file.  Prints resulting gff3 to STDOUT.

=head1 USAGE

    gff3_insert_sync_marks.pl  file.gff3  anotherfile.gff3 ... > with_syncs.gff3

=cut
