#!/bin/bash

perl -e '
    use strict;
    use warnings;

    use LinearAlgebra qw( mult_matrix_product );

    my $matrix1 = [ [ 1 ], [ 0 ], [ 0 ], ];
    my $matrix2 = [ [ 3 ], [ 2 ], [ 1 ], ];

    my $matrix_product = mult_matrix_product( [ $matrix1, $matrix2 ] );

    for my $matrix_id ( 0..$#{ $matrix_product } ) {
        for my $row ( @{ $matrix_product->[$matrix_id] } ) {
            print( join( " ", @{ $row } ), "\n" );
        }
        print( "\n" ) if $matrix_id != $#{ $matrix_product };
    }
' 2>&1 | sed 's/line\s*[0-9]*.$/line <row>./g' | sed 's/0x[0-9a-f]\{12\}/<hex>/g'
