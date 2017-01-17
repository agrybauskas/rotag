package LinearAlgebra;

use strict;
use warnings;

use Math::Algebra::Symbols;

# Describes path to Maxima interpreter tool.
# TODO: look, how to do declare path variable properly.
my $maxima_path = "../../tools";

# ------------------------------- Linear algebra ------------------------------ #

#
# Constants
#

my $PI = 4 * atan2( 1, 1 );
my $EPSILON = 1.0; # Machine accuracy for floating point numbers which is
                   # calculated in above block. 

{
    while( ( 1.0 + 0.5 * $EPSILON ) != 1.0 ) {
	$EPSILON = 0.5 * $EPSILON;
    }
}

#
# Creates local reference frame for any three given atoms positions in cartesian
# coordinate system.
# Input  (1 arg): array of three atom coordinates in x, y, z form.
# Output (1 arg): array of reference frame coordinates in x, y, z form.
#

sub create_ref_frame
{
    my ( $mid_atom_x,  $mid_atom_y,  $mid_atom_z,
         $up_atom_x,   $up_atom_y,   $up_atom_z,
         $side_atom_x, $side_atom_y, $side_atom_z ) = @_;

    my @local_ref_frame;

    # Let local z-axis be colinear to bond between mid and up atoms.
    $local_ref_frame[2][0] = $up_atom_x - $mid_atom_x;
    $local_ref_frame[2][1] = $up_atom_y - $mid_atom_y;
    $local_ref_frame[2][2] = $up_atom_z - $mid_atom_z;

    # Let local x-axis be perpendicular to bonds between mid, up and mid, side
    # atoms.
    $local_ref_frame[0][0] =
        ( $side_atom_y - $mid_atom_y ) * $local_ref_frame[2][2]
      - ( $side_atom_z - $mid_atom_z ) * $local_ref_frame[2][1];
    $local_ref_frame[0][1] =
      - ( $side_atom_x - $mid_atom_x ) * $local_ref_frame[2][2]
      + ( $side_atom_z - $mid_atom_z ) * $local_ref_frame[2][0];
    $local_ref_frame[0][2] =
        ( $side_atom_x - $mid_atom_x ) * $local_ref_frame[2][1]
      - ( $side_atom_y - $mid_atom_y ) * $local_ref_frame[2][0];

    # Let local y-axis be in the same plane as mid-up and mid-side bonds.
    $local_ref_frame[1][0] =
        $local_ref_frame[2][1] * $local_ref_frame[0][2]
      - $local_ref_frame[2][2] * $local_ref_frame[0][1];
    $local_ref_frame[1][1] =
      - $local_ref_frame[2][0] * $local_ref_frame[0][2]
      + $local_ref_frame[2][2] * $local_ref_frame[0][0];
    $local_ref_frame[1][2] =
        $local_ref_frame[2][0] * $local_ref_frame[0][1]
      - $local_ref_frame[2][1] * $local_ref_frame[0][0];

    return @local_ref_frame;
}

#
# Function calculates Euler rotational angles (alpha, beta, gamma) that are used
# to transform global reference frame to chosen one.
# Input  (1 arg): array of three atom coordinates in x, y, z form.
# Output (3 arg): euler angles (alpha, beta, gamma) in radians.
#

sub find_euler_angles
{
    my ( $mid_atom_x,  $mid_atom_y,  $mid_atom_z,
         $up_atom_x,   $up_atom_y,   $up_atom_z,
         $side_atom_x, $side_atom_y, $side_atom_z ) = @_;

    my $alpha_rad;
    my $beta_rad;
    my $gamma_rad;

    my $z_axis_in_xy_plane;

    my @local_ref_frame =
        create_ref_frame( $mid_atom_x,  $mid_atom_y,   $mid_atom_z,
                          $up_atom_x,   $up_atom_y,    $up_atom_z,
                          $side_atom_x, $side_atom_y,  $side_atom_z );

    # Projects local z-axis to global xy-plane.
    $z_axis_in_xy_plane =
        sqrt( $local_ref_frame[2][0] * $local_ref_frame[2][0]
            + $local_ref_frame[2][1] * $local_ref_frame[2][1] );

    if( $z_axis_in_xy_plane > $EPSILON ) {
        $alpha_rad =
            atan2( $local_ref_frame[1][0] * $local_ref_frame[2][1]
                 - $local_ref_frame[1][1] * $local_ref_frame[2][0],
                   $local_ref_frame[0][0] * $local_ref_frame[2][1]
                 - $local_ref_frame[0][1] * $local_ref_frame[2][0] );
        $beta_rad = atan2( $z_axis_in_xy_plane, $local_ref_frame[2][2] );
        $gamma_rad = - atan2( - $local_ref_frame[2][0], $local_ref_frame[2][1] );
    } else {
        $alpha_rad = 0.;
        $beta_rad = ( $local_ref_frame[2][2] > 0. ) ? 0. : $PI;
        $gamma_rad = - atan2( $local_ref_frame[0][1], $local_ref_frame[0][0] );
    }

    return $alpha_rad, $beta_rad, $gamma_rad;
}

#
# Transposes matrix.
# Input:  (1 arg): array representing matrix.
# Output: (1 arg): transposed matrix.
#

sub transpose
{
    my $matrix = shift;
    my @matrix = @$matrix;

    my @transposed_matrix;

    for my $row ( 0..$#matrix ) {
	for my $col ( 0..$#{ $matrix[$row] } ) {
	    $transposed_matrix[$col][$row] = $matrix[$row][$col];
	}
    }

    return \@transposed_matrix;
}

# ---------------------------- Symbolic linear algebra ------------------------ #

#
# Example of rotation along z-axis by chi angle in radians:
#
#      / cos(chi) -sin(chi) 0 \   / x \   / x * cos(chi) + y * sin(chi) \
#      | sin(chi)  cos(chi) 0 | * | y | = | x * sin(chi) + y * cos(chi) |
#      \    0         0     1 /   \ z /   \              z              /
#

#
# Calculates dot product of two matrices that might have symbolic variables.
# Input:  (3 arg): first argument - symbols that identify strings as
#                  symbols for mathematical manipulation, second and third - 2
#                  arrays that are correctly paired.
# Output: (1 arg): dot product.
#

sub two_matrix_product
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @matrix_product;

    # Notifies error, when the column number of left matrix does not equal the
    # row number of the right matrix.
    die( "A row number of a left matrix is NOT equal to the column\n" .
    	 "number of the right matrix.\n" )
    	unless( scalar( @{ transpose( $left_matrix ) } ) ==
    		scalar( @$right_matrix ) );

    # Makes placeholder items for matrix_product array.
    for( my $product_row = 0; 
    	 $product_row < scalar( @$left_matrix ); 
    	 $product_row++ ) {
    	for( my $product_col = 0; 
    	     $product_col < scalar( @{ $right_matrix->[0] } ); 
    	     $product_col++ ) {
    	    $matrix_product[$product_row][$product_col] = 0;
    	}
    }

    # Calculates dot product.
    for( my $product_row = 0;
    	 $product_row < scalar( @matrix_product );
    	 $product_row++ ) {
    	for( my $product_col = 0;
    	     $product_col < scalar( @{ $matrix_product[$product_row] } );
    	     $product_col++ ) {
    	    for( my $left_col = 0;
    		 $left_col < scalar( @{ $left_matrix->[$product_col] } );
    		 $left_col++ ) {
    		my $left_number;
    		my $right_number;
		my $dot_command;

		$left_number = $left_matrix->[$product_row]->[$left_col];
		$right_number = $right_matrix->[$left_col]->[$product_col];

		$dot_command = 
		      "$maxima_path/maxima_listener \"fortran(("
		      . "$matrix_product[$product_row][$product_col]"
		      . ")+($left_number)*($right_number));\"";

		$matrix_product[$product_row][$product_col] = qx( $dot_command );
	    }
    	}
    }

    return \@matrix_product;
}

#
# TODO: write full discription about a function.
#

sub mult_matrix_product
{
    my @matrices = @_;

    my $mult_matrix_product;

    for( my $id = $#matrices; $id >= 1; $id-- ) {
	if( $id == $#matrices ) {
    	    $mult_matrix_product = two_matrix_product( $matrices[$id-1],
						       $matrices[$id] );
    	} else {
    	    $mult_matrix_product = two_matrix_product( $matrices[$id-1],
						       $mult_matrix_product );
    	}
    }

    return $mult_matrix_product;
}

1;
