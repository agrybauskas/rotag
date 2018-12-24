#!/usr/bin/perl

use strict;
use warnings;

use LinearAlgebra qw( mult_matrix_product );
use Symbolic;

my $matrix1 = [ [ 1, 0, 0 ],
                [ 0, 2, 0 ],
                [ 0, 0, 3 ], ];
my $matrix2 =
    Symbolic->new(
        { 'symbols' => [ 'x', 'y', 'z' ],
          'matrix' => sub { my ( $x, $y, $z ) = @_;
                            return [ [ $x ],
                                     [ $y ],
                                     [ $z ], ] } } );

my $matrix_product =
    mult_matrix_product( [ $matrix1, $matrix2 ],
                         { 'x' => 6, 'y' => 3, 'z' => 2 } );

for my $matrix_id ( 0..$#{ $matrix_product } ) {
    for my $row ( @{ $matrix_product->[$matrix_id] } ) {
        print( join( " ", @{ $row } ), "\n" );
    }
    print( "\n" ) if $matrix_id != $#{ $matrix_product };
}
