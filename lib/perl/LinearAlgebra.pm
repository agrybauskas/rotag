package LinearAlgebra;

use Math::Algebra::Symbols;

use strict;
use warnings;

# ------------------------------- Linear algebra ------------------------------ #

#
# Block of code contains functions that perform basic linear algebra operations
# for matrices.
#

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
# Switches between global and local reference frames.
# Input
#

sub switch_ref_frame
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $switch_to_local ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	LinearAlgebra::find_euler_angles( @$mid_atom_coord,
					  @$up_atom_coord,
					  @$side_atom_coord );

    # Depending on the option switch_to_local, 
    my @switch_matrix;

    if( $switch_to_local == 1 ) {
	@switch_matrix =
	    LinearAlgebra::mult_matrix_product(
		LinearAlgebra::rotate_z_axis( $alpha ),
		LinearAlgebra::rotate_x_axis( $beta ),
		LinearAlgebra::rotate_z_axis( $gamma ),
		LinearAlgebra::translate( ( - $mid_atom_coord->[0],
					    - $mid_atom_coord->[1],
					    - $mid_atom_coord->[2] ) ) )
    } elsif( $switch_to_local == 0 ) {
	@switch_matrix =
	    LinearAlgebra::mult_matrix_product(
		LinearAlgebra::rotate_z_axis( - $alpha ),
		LinearAlgebra::rotate_x_axis( - $beta ),
		LinearAlgebra::rotate_z_axis( - $gamma ),
		LinearAlgebra::translate( ( $mid_atom_coord->[0],
					    $mid_atom_coord->[1],
					    $mid_atom_coord->[2] ) ) )
    } else {
	die "Must choose \$switch_to_global value between 0 and 1.\n"
    }

    return \@switch_matrix;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around x-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub rotate_x_axis
{
    my $angle = shift;

    my @rot_matrix_x =
	( [ 1, 0, 0, 0 ],
	  [ 0, cos( $angle ), - sin( $angle ), 0 ],
	  [ 0, sin( $angle ),   cos( $angle ), 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_x;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around y-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub rotate_y_axis
{
    my $angle = shift;

    my @rot_matrix_y =
	( [   cos( $angle ), 0, sin( $angle ), 0 ],
	  [ 0, 1, 0, 0 ],
	  [ - sin( $angle ), 0, cos( $angle ), 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_y;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around z-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub rotate_z_axis
{
    my $angle = shift;

    my @rot_matrix_z =
	( [ cos( $angle ), - sin( $angle ), 0, 0 ],
	  [ sin( $angle ),   cos( $angle ), 0, 0 ],
	  [ 0, 0, 1, 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_z;
}

#
# Creates 4x4 matrix that can translates 4x4 matrices
# Input  (1 arg): 3x1 matrix of x, y, z coordinates that .
# Output (1 arg): 4x4 matrix.
#

sub translate
{
    my @transl_coord = @_;

    my @transl_matrix =
	( [ 1, 0, 0, $transl_coord[0] ],
	  [ 0, 1, 0, $transl_coord[1] ],
	  [ 0, 0, 1, $transl_coord[2] ],
	  [ 0, 0, 0, 1 ] );

    return \@transl_matrix;
}

# ---------------------------- Symbolic linear algebra ------------------------ #

#
# Block of code contains functions that perform basic linear algebra on symbolic
# expressions..
#
#
# Example of rotation along z-axis by chi angle:
#
#      / cos(chi) -sin(chi) 0 \   / x \   / x * cos(chi) + y * sin(chi) \
#      | sin(chi)  cos(chi) 0 | * | y | = | x * sin(chi) + y * cos(chi) |
#      \    0         0     1 /   \ z /   \              z              /
#

#
# Transposes matrix.
# Input  (1 arg): array representing matrix.
# Output (1 arg): transposed matrix.
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

#
# Calculates matrix product of two matrices that might have symbolic variables.
# Input  (2 arg): 2 arrays that are correctly paired.
# Output (1 arg): matrix product.
#

sub two_matrix_product
{
    my ( $symbols, $left_matrix, $right_matrix ) = @_;

    my %symbols; # Hash that prepares symbols for algebraic manipulation
    my @matrix_product;

    # Notifies error, when the column number of left matrix does not equal the
    # row number of the right matrix.
    die(   "A row number of a left matrix is NOT equal to the column\n"
         . "number of the right matrix.\n" )
    	unless( scalar( @{ transpose( $left_matrix ) } ) ==
    		scalar( @$right_matrix ) );

    # Makes placeholder items consisting zero values for matrix_product array.
    for( my $product_row = 0;
    	 $product_row < scalar( @$left_matrix );
    	 $product_row++ ) {
    	for( my $product_col = 0;
    	     $product_col < scalar( @{ $right_matrix->[0] } );
    	     $product_col++ ) {
    	    $matrix_product[$product_row][$product_col] = 0;
    	}
    }

    # Initiates perception of symbols.
    foreach( @$symbols ) {
	$symbols{$_} = symbols( $_ );
    }

    # Calculates matrix product of two matrices.
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

		# Retrieves numbers that will be multiplied and added to
		# matrix_product array.
		$left_number = $left_matrix->[$product_row][$left_col];
		$right_number = $right_matrix->[$left_col][$product_col];

		# Changes "$" to hash reference (for symbols that are
		# written in "$x" form).
		$left_number =~ s/\$(\w+)/\$symbols{$1}/g;
		$right_number =~ s/\$(\w+)/\$symbols{$1}/g;

		# Changes "&" to hash reference (for symbols that are
		# written in "&i" form, usually for Euler number).
		$left_number =~ s/\&(\w+)/\$symbols{$1}/g;
		$right_number =~ s/\&(\w+)/\$symbols{$1}/g;

		$matrix_product[$product_row][$product_col] +=
		    eval( $left_number ) * eval( $right_number );
	    }
    	}
    }

    return \@matrix_product;
}

#
# Calculates matrix product of list of any size of matrices.
# Input  (n arg): any number of arrays representing matrices.
# Output (1 arg): matrix product.
#

sub mult_matrix_product
{
    my $symbols = shift;
    my @matrices = @_;

    my $mult_matrix_product;

    # Multiplies matrices from left to right.
    for( my $id = $#matrices; $id >= 1; $id-- ) {
	if( $id == $#matrices ) {
    	    $mult_matrix_product = two_matrix_product( $symbols,
						       $matrices[$id-1],
						       $matrices[$id] );
    	} else {
    	    $mult_matrix_product = two_matrix_product( $symbols,
						       $matrices[$id-1],
						       $mult_matrix_product );
    	}
    }

    return $mult_matrix_product;
}

1;
