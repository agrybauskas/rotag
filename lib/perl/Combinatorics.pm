package Combinatorics;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( permutation );
use Data::Dumper;
# ---------------------------------- Combinatorics ---------------------------- #

#
# Carries out permutation for list of arrays.
#

# TODO: move global variable inside function.
my @permuted_list;

sub permutation
{
    my ( $size, $base, $list ) = @_;
    my $base_size = scalar( @{ $base } );

    if( $base_size == $size ) {
	push( @permuted_list, $base );
    } else {
	for my $i ( @{ $list->[$base_size] } ) {
	    permutation( $size, [ @{ $base }, $i ], $list );
	}
    }

    return \@permuted_list;
}

1;
