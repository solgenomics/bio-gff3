use strict;
use warnings;

use File::Temp;
use IO::String;
use Test::More;

use Bio::GFF3::Transform::FromFasta 'gff3_from_fasta';

my $test1_fasta = <<'';
>foo  and this is my description  
ACTGA  TGATCTGATGATGATGCTTAG TGCTGA
GATCTGATGAC
> conan  The Barbarian!!
CTAAGATGCTCGATGATAGCTAGCATGATGCTAGCTAGCATGCTAGCATGATGCATCGATCATGATG

my $test1_gff3 = _platform_lines(<<'');
##gff-version 3
foo	fasta	region	1	44	.	+	.	Name=foo;Note=and this is my description  
conan	fasta	region	1	67	.	+	.	Name=conan;Note=The Barbarian!!
##FASTA
>foo  and this is my description  
ACTGATGATCTGATGATGATGCTTAGTGCTGA
GATCTGATGAC
> conan  The Barbarian!!
CTAAGATGCTCGATGATAGCTAGCATGATGCTAGCTAGCATGCTAGCATGATGCATCGATCATGATG

my $out;
gff3_from_fasta(
    in   => \$test1_fasta,
    out  => \$out,
    type => 'region',
    fasta_section => 1,
  );

is $out, $test1_gff3, 'right out';

my $tempfile = File::Temp->new;
$tempfile->print( $test1_fasta );
$tempfile->close;

my $out2;
gff3_from_fasta(
    in   => "$tempfile",
    out  => \$out2,
    type => 'region',
    fasta_section => 1,
  );


is $out, $test1_gff3, 'right out reading from tempfile';

done_testing;

# add windows line endings to test strings if necessary
sub _platform_lines {
    if( $^O eq 'MSWin32' ) {
        s/\n/\r\n/g for @_;
    }
    return wantarray ? @_  : $_[0]
}
