package Combinatorics;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( permutation );

# ---------------------------------- Combinatorics ---------------------------- #

#
# Carries out permutation for list of arrays (recursive function).
# Input:
#     $size - size of a list.
#     $base - list of items that permutation list will be generated from.
#     $list - one item of whole permutation product.
#     $permuted_list - permutation product that is tracked through all
#     recursions.
# Output:
#     $permuted_list - final permutation product.
#

sub permutation
{
    my ( $size, $base, $list, $permuted_list ) = @_;
    my $base_size = scalar( @{ $base } );

    if( $base_size == $size ) {
	push( @{ $permuted_list }, $base );
    } else {
	for my $i ( @{ $list->[$base_size] } ) {
	    permutation( $size, [ @{ $base }, $i ], $list, $permuted_list );
	}
    }

    return $permuted_list;
}

1;
