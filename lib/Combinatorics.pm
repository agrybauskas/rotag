package Combinatorics;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( adjacent
                     permutation );

use List::Util qw( max );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ------------------------------ Combinatorics -------------------------------- #

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
    my $base_size = scalar @{ $base };

    if( $base_size == $size ) {
        push @{ $permuted_list }, $base;
    } else {
        for my $i ( @{ $list->[$base_size] } ) {
            permutation( $size, [ @{ $base }, $i ], $list, $permuted_list );
        }
    }

    return $permuted_list;
}

sub adjacent
{
    my ( $lists ) = @_;

    # Find list with max items.
    my $max_item_count = max( map { scalar( @{ $_ } ) } @{ $lists } );

    my @list = ();
    for my $i ( 0..$max_item_count-1 ) {
        my @current_list = ();
        # If the list has less items than the list that has the most items,
        # it circulary chooses next.
        for my $list ( @{ $lists } ) {
            my $item_count = scalar @{ $list };
            my $i_updated = ( ( $i + 1 ) % $item_count ) - 1;
            push @current_list, $list->[$i_updated];
        }
        push @list, \@current_list;
    }

    return \@list;
}

1;
