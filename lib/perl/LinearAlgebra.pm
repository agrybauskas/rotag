package LinearAlgebra;

use strict;
use warnings;

use Data::Dumper;
use feature qw( say );

use Math::Algebra::Symbols;

# ------------------------------- Linear algebra ------------------------------ #

#
# Constants
#

my $PI = 4 * atan2( 1, 1 );
my $EPSILON = 1.0 / ( 2 ** 52 ); # Machine accuracy for 64-bit floating point
                                 # numbers.
                                 # TODO: make machine accuracy dependent on the
                                 # machine.

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
# Calculates dot product of two matrices.
# Input:  (2 arg): 2 arrays that are correctly paired.
# Output: (1 arg): dot product.
#

sub dot_product
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @dot_product;

    # TODO: make error notification, when the column number of left matrix
    #       does not equal row number of right matrix.

    # Makes placeholder items for dot_product array.
    for( my $dot_row = 0; 
	 $dot_row < scalar( @$left_matrix ); 
	 $dot_row++ ) {
	for( my $dot_col = 0; 
	     $dot_col < scalar( @{ $right_matrix->[0] } ); 
	     $dot_col++ ) {
	    $dot_product[$dot_row][$dot_col] = 0;
	}
    }

    # Calculates dot product.
    for( my $dot_row = 0;
	 $dot_row < scalar( @dot_product );
	 $dot_row++ ) {
	for( my $dot_col = 0;
	     $dot_col < scalar( @{ $dot_product[$dot_row] } );
	     $dot_col++ ) {
	    for( my $left_col = 0;
		 $left_col < scalar( @{ $left_matrix->[$dot_col] } );
		 $left_col++ ) {
		$dot_product[$dot_row][$dot_col] += 
		    eval( $left_matrix->[$dot_row]->[$left_col] ) 
		  * eval( $right_matrix->[$left_col]->[$dot_col] );
	    }
	}
    }

    return \@dot_product;
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
#      \    0         0     1 /   \ z /   \              0              /
#

#
# Calculates dot product of given matrices recursively from right to left. All
# matrices has to be in correct size in order to be multiplied correctly.
# Input  (2 to n arg): first argument - symbols that identify strings as
#                      symbols for mathematical manipulation, others - any
#                      number of arrays of arrays representing n x n matrices. 
# Output      (1 arg): array representing correctly calculated dot product.
#

sub symb_dot_product
{
    my ( $symbols, $symb_matrices ) = @_;

    my @symb_matrices = @$symb_matrices;
    my %symbols; # Hash that prepares symbols for algebraic manipulation.

    foreach( @$symbols ) {
    	$symbols{$_} = symbols( $_ );
    }
    
}

1;
