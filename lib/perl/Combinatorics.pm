package Combinatorics;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( permutation );

# ---------------------------------- Combinatorics ---------------------------- #

#
# Carries out permutation for list of arrays.
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
