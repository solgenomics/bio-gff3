package Bio::GFF3::Transform::Reformat;
use strict;
use warnings;
use Carp;

use Data::Dumper;

use Tie::Cache::LRU;
use URI::Escape;

my $default_uniq_sep = '_';
my $default_lb_size = 200;

sub reformat {
    my ( $out_fh, %args ) = @_;

    tie my %uniq_lb, 'Tie::Cache::LRU', $lb_size;
    my %uniq_ctrs;

    #if we have a -f option, read the filenames out of that file and add them to ARGV
    if( $opt{f} ) {
        open my $files, $opt{f} or die "$! reading $opt{f}";
        while( <$files> ) {
            s/^\s+|\s+$//g;
            -f && -r or die "file $_ is not readable\n";
            push @ARGV,$_;
        }
    }


    my $seqregions = do {
        if($opt{S}) {
            -r $opt{S} or die "cannot open '$opt{S}' for reading\n";
            index_seqlengths($opt{S});
        } else {
            {}
        }
    };

    print $out_fh "##gff-version 3\n";

    unless ( $opt{i} ) { #< if we are not interleaving seqregions, print them all here
        foreach my $sr_rec (map $seqregions->{$_}, sort keys %$seqregions) {
            print $out_fh sr_str($sr_rec)."\n";
            $sr_rec->{printed} = 1;
        }
    }

    my $in_fh;
    if ( $opt{s} ) {
        open $in_fh, "sort -k 1,1 -k 4,4g -s @ARGV | grep -v '^###' |"
            or die "$! running sort";
    } else {
        open $in_fh, "cat @ARGV |"
            or die "$! running cat";
    }

    while( my $line = <$in_fh> ) {
        $line = $opt{L}->($line) if $opt{L}; #< do the global -L alteration if present
        next if $line =~ /^##gff-version/;

        if ( $line =~ /^##\s*sequence-region/ ) {

            #check sequence-region directive, check but don't repeat
            #directives that have already been printed
            chomp $line;
            my (undef,$seqname,$start,$end) = split /\s+/,$line;

            if ( my $known_sr = $seqregions->{$seqname} ) { #< if we already know about this sequence-region, check it and print it if 
                my $this_sr = { name => $seqname, start => $start, end => $end, length => $end-$start+1 };
                sr_eq( $this_sr, $known_sr )
                    or warn "WARNING: sequence-region statement '".sr_str($this_sr)."' conflicts with previously-seen sequence region length (".sr_str($known_sr).") from -S file or earlier in the GFF stream.  Overriding with first-seen length.\n";

            } else {
                $seqregions->{$seqname} = {
                    length  => $end-$start+1,
                    start   => $start,
                    end     => $end,
                    name    => $seqname,
                    printed => 0,
                };
            }
        } elsif ( $line =~ /^\S+\t\S+\t\S+/ ) { #a data line, process it
            chomp $line;
            my @fields = split /\s+/, $line, 9;
            my $fcnt = @fields;
            $fcnt == 9 or die "invalid number of fields ($fcnt)";

            if ( my $sr_rec = $seqregions->{$fields[0]} ) {
                unless ( $sr_rec->{printed} ) {
                    print $out_fh sr_str($sr_rec)."\n";
                    $sr_rec->{printed} = 1;
                }
            }

            #use Data::Dumper
            #warn Dumper \@attrs;

            # if we have -U or attr expressions, we have to parse and mess with the attributes
            if ( $opt{U} || @attr_exprs ) {

                #parse the attributes
                my @attrs = map [split /=/],split /;/,$fields[8];

                foreach my $a (@attrs) {
                    my ($name,$val) = @$a;
                    $val = uri_unescape($val);

                    foreach my $ae (@attr_exprs) {
                        my ($qr,$change_sub) = @$ae;
                        #warn "try matching $name with '$qr' => '$expr'\n";
                        if ( $name =~ $qr) {
                            #warn "matched $qr, $change_sub\n";
                            $val = $change_sub->($val);
                        }
                    }

                    if ( $opt{U} ) {
                        if ( $name eq 'Parent') {
                            my $key = join ':',@fields[0,1],$val;
                            #find the uniqified version in the lookback buffer
                            $uniq_lb{$key} or die "no feature found with key '$key', either this file is not valid GFF3, or you have parent and child features very far away from eachother in this file and need to increase the lookback buffer size with the -l option (currently -l $lb_size).\nCurrent lookback buffer contents: ".Dumper(\%uniq_lb);
                            $val = $uniq_lb{$key};
                        } elsif ( $name eq 'ID' ) {
                            my $new = $val;
                            $new =~ s/$uniq_sep\d+$//;
                            my $index = ++$uniq_ctrs{$new};
                            $new .= $uniq_sep.$index unless $index == 1;

                            my $key = join ':', @fields[0,1], $val;

                            unless ( $opt{e} && $new =~ $opt{e} ) {
                                $uniq_lb{$key} = $new;
                            }
                            $val = $new;
                        }
                    }

                    # add an ID if we got -I option and this feature has no ID
                    if ( $opt{I} ) {
                        unless( grep $_->[0] eq 'ID', @attrs ) {
                            unshift @attrs,['ID',$fields[0].'_'.$fields[2].'_'.++$uniq_ctrs{$fields[2]}];
                        }

                        # note: we don't have to worry about Parent attrs, because
                        # anything that doesn't have an ID will not have any
                        # elements referring to it as Parent
                    }

                    $a = [$name,uri_escape($val,"\t\n".';=%&,[:cntrl:]')]; #<alter the attributes
                }

                $fields[8] = join ';', map join('=',@$_), @attrs;
            }

            print $out_fh join("\t",@fields);
            print $out_fh "\n";
        } else { 		#some other thing, just print it
            print $out_fh $_;
        }
    }
}

#given a sequence file, return a hashref of its sequence lengths
sub index_seqlengths {
  my $seqfile = shift;

  require Bio::SeqIO;
  my $seq_in = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');
  my %lengths;
  while( my $s = $seq_in->next_seq ) {
    $lengths{$s->primary_id} = { start => 1, end => $s->length, length => $s->length, name => $s->primary_id, printed => 0 }
      or die "in '$seqfile', sequence ".$s->primary_id." has no length!\n";
  }
  return \%lengths;
}

# couple of small functions for dealing with sequence-region records
sub sr_eq {
  my ($one,$two) = @_;
  foreach (qw/length start end/) {
    return 0 unless $one->{$_} == $two->{$_};
  }
  return 1;
}
sub sr_str {
  my ($sr) = @_;
  return "##sequence-region $sr->{name} $sr->{start} $sr->{end}";
}

1;
