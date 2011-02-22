use strict;
use warnings;

use Test::More;

use Bio::GFF3::LowLevel qw(
    gff3_parse_feature
    gff3_parse_attributes
    gff3_format_feature
    gff3_format_attributes
    gff3_escape
    gff3_unescape
);

roundtrip(
    "FooSeq\tbarsource\tmatch\t234\t234\t0.0\t+\t.\tID=Beep%2Cbonk%3B+Foo\n",
    {
        'attributes' => {
            'ID' => [
                'Beep,bonk;+Foo'
              ]
          },
        'end' => '234',
        'phase' => undef,
        'score' => '0.0',
        'seqid' => 'FooSeq',
        'source' => 'barsource',
        'start' => '234',
        'strand' => '+',
        'type' => 'match'
      },
  );

roundtrip(
    gff3_escape("Noggin,+-\%Foo\tbar")."\tbarsource\tmatch\t234\t234\t0.0\t+\t.\t.\n",
    {
        'attributes' => undef,
        'end' => '234',
        'phase' => undef,
        'score' => '0.0',
        'seqid' => "Noggin,+-\%Foo\tbar",
        'source' => 'barsource',
        'start' => '234',
        'strand' => '+',
        'type' => 'match'
    },
  );

done_testing;


####

sub roundtrip {
    my ( $line, $parsed ) = @_;

    my $p = gff3_parse_feature( $line );
    is_deeply( $p, $parsed, "parsed correctly" )
      or diag explain $p;

    is( gff3_format_feature($parsed), $line, 'round-tripped back to original line' );

}
