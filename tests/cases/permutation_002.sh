#!/usr/bin/perl

use strict;
use warnings;

use Combinatorics qw( permutation );

my $permuted_list = permutation( 2,
                                 [],
                                 [ [ [ 0 ], [ 1 ], [ 2 ] ],
                                   [ [ 3, 5 ], [ 6, 7 ] ] ],
                                 [] );
$permuted_list = [ map { [ @{ $_->[0] }, @{ $_->[1] } ] } @{ $permuted_list } ];

print( join( " ", @{ $_ } ), "\n" ) for @{ $permuted_list };
