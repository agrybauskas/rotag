package LinearAlgebra;

use Exporter qw( import );
@EXPORT_OK = qw( create_ref_frame
                 find_euler_angles
                 switch_ref_frame
                 x_axis_rotation
                 y_axis_rotation
                 z_axis_rotation
                 translation
                 vector_cross
                 matrix_sum
                 matrix_sub
                 transpose
                 matrix_product
                 evaluate_matrix );

use strict;
use warnings;
use Data::Dumper;
# --------------------------- Numeric linear algebra -------------------------- #

#
# Performs basic linear algebra operations for matrices.
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

    return @local_ref_frame; # TODO: return to reference.
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
# Switches between global and local frames of reference.
# Input  (2 arg): array of three atom coordinates in x, y, z form, that frame of
#                 reference will be transformed to, and boolean that switches
#                 frame of reference to local (if "local") and global
#                 (if "global").
# Output (1 arg): 4x4 transformation matrix.
#

sub switch_ref_frame
{
    my ( $switch_to_local,
	 $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) = find_euler_angles( @$mid_atom_coord,
						      @$up_atom_coord,
						      @$side_atom_coord );
    # Depending on the option switch_to_local,
    my @switch_matrix;

    if( $switch_to_local eq "local" ) {
    	@switch_matrix =
    	    matrix_product( &z_axis_rotation( $alpha ),
    			    &x_axis_rotation( $beta ),
    			    &z_axis_rotation( $gamma ),
    			    &translation( ( -$mid_atom_coord->[0],
					    -$mid_atom_coord->[1],
					    -$mid_atom_coord->[2] ) ) );
    } elsif( $switch_to_local eq "global" ) {
    	@switch_matrix =
    	    matrix_product( &translation( ( $mid_atom_coord->[0],
					    $mid_atom_coord->[1],
					    $mid_atom_coord->[2] ) ),
    			    &z_axis_rotation( -$gamma ),
    			    &x_axis_rotation( -$beta ),
    			    &z_axis_rotation( -$alpha ) );
    } else {
    	die "Must choose \$switch_to_global value between \"local\" and"
    	    . "\"global\".\n"
    }

    return @switch_matrix;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around x-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub x_axis_rotation
{
    my ( $angle ) = @_;

    my @rot_matrix_x =
	( [ 1, 0, 0, 0 ],
	  [ 0, cos( $angle ), -sin( $angle ), 0 ],
	  [ 0, sin( $angle ),  cos( $angle ), 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_x;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around y-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub y_axis_rotation
{
    my ( $angle ) = @_;

    my @rot_matrix_y =
	( [  cos( $angle ), 0, sin( $angle ), 0 ],
	  [ 0, 1, 0, 0 ],
	  [ -sin( $angle ), 0, cos( $angle ), 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_y;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around z-axis by making a
# matrix product.
# Input  (1 arg): angle of rotation in radians.
# Output (1 arg): 4x4 matrix.
#

sub z_axis_rotation
{
    my ( $angle ) = @_;

    my @rot_matrix_z =
	( [ cos( $angle ), -sin( $angle ), 0, 0 ],
	  [ sin( $angle ),  cos( $angle ), 0, 0 ],
	  [ 0, 0, 1, 0 ],
	  [ 0, 0, 0, 1 ] );

    return \@rot_matrix_z;
}

#
# Creates 4x4 matrix that can translates 4x4 matrices by making a matrix product.
# Input  (1 arg): 3x1 matrix of x, y, z coordinates of corresponding
#                 displacement.
# Output (1 arg): 4x4 matrix.
#

sub translation
{
    my @transl_coord = @_;

    my @transl_matrix =
	( [ 1, 0, 0, $transl_coord[0] ],
	  [ 0, 1, 0, $transl_coord[1] ],
	  [ 0, 0, 1, $transl_coord[2] ],
	  [ 0, 0, 0, 1 ] );

    return \@transl_matrix;
}

#
# Produces cross product of two 3D vectors.
# Input  (2 arg): two arrays representing 3D vectors where units are i, j, k.
# Output (1 arg): cross product.
#

sub vector_cross
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @cross_product =
	( $left_matrix->[1] * $right_matrix->[2]
	 -$left_matrix->[2] * $right_matrix->[1],
	 -$left_matrix->[0] * $right_matrix->[2]
	 +$left_matrix->[2] * $right_matrix->[0],
	  $left_matrix->[0] * $right_matrix->[1]
	 -$left_matrix->[1] * $right_matrix->[0] );

    return \@cross_product;
}

#
# Adds two matrices.
# Input  (2 arg): two matrices.
# Output (1 arg): sum of matrices.
#

sub matrix_sum
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @matrix_sum;

    for my $i ( 0..$#{ $left_matrix } ) {
    	for my $j ( 0..$#{ $left_matrix->[$i] } ) {
	    $matrix_sum[$i][$j] =
		$left_matrix->[$i][$j] + $right_matrix->[$i][$j];
	}
    }

    return \@matrix_sum;
}

#
# Substracts one matrix from another.
# Input  (2 arg): two matrices.
# Output (1 arg): difference of matrices.
#

sub matrix_sub
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @matrix_sub;

    for my $i ( 0..$#{ $left_matrix } ) {
    	for my $j ( 0..$#{ $left_matrix->[$i] } ) {
	    $matrix_sub[$i][$j] =
		$left_matrix->[$i][$j] - $right_matrix->[$i][$j];
	}
    }

    return \@matrix_sub;
}

# ---------------------------- Symbolic linear algebra ------------------------ #

#
# Performs basic linear algebra on symbolic expressions. Uses GiNaC for
# calculations.
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
    my ( $matrix ) = @_;
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
# Calculates matrix product of list of any size of matrices.
# Input  (n arg): any number of arrays representing matrices.
# Output (1 arg): matrix product.
#

sub matrix_product
{
    my @matrices = @_;

    # Converts matrices to GiNaC readable input.
    my $matrix_equation_ginac =
	join( "*",
	map { "[" . join( ",",
	map { "[" . join( ",",
	@$_ ) . "]" }
        @$_ ) . "]" }
	@matrices );

    # Converts perl power symbol ** to GiNaC's ^.
    $matrix_equation_ginac =~ s/\*\*/^/g;

    # Runs GiNaC.
    my $matrix_product_ginac =
	qx/ echo "evalm( ${matrix_equation_ginac} );" | ginsh /
	|| die( "A row number of a left matrix is NOT equal to the column\n" .
		"number of the right matrix.\n" );

    # Returns matrices in perl array form.
    my @matrix_product;

    for my $row ( split( ",\\[", $matrix_product_ginac ) ) {
	$row =~ s/\n//g; # Remove newline.
	$row =~ s/]//g; # Removes unnecessary GiNaC matrix symbols.
	$row =~ s/\[//g; # Removes unnecessary GiNaC matrix symbols.
	$row =~ s/\^/\*\*/g; # Brings back perl power symbol **.
	push( @matrix_product, [ split( ",", $row ) ] );
    }

    return \@matrix_product;
}

#
# Evaluates symbolic variables.
# Input  (1 arg): matrix with variables as symbols.
# Output (1 arg): evaluated matrix.
#

sub evaluate_matrix {
    my ( $symbols, $matrix_ginac ) = @_;
    my %symbols = %$symbols;

    my @eval_matrix;

    # Adds $ sign to given symbols and then runs eval() function.
    for my $row ( @$matrix_ginac ) {
	push( @eval_matrix, [] );
	for my $value ( @$row ) {
	    for my $symbol ( keys %$symbols ) {
		$value =~ s/\b(${symbol})\b/\$symbols{$1}/g;
	    }
	    push( @{ $eval_matrix[-1] }, eval( $value ) );
	}
    }

    return \@eval_matrix;
}

1;
