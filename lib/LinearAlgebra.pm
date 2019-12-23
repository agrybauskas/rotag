package LinearAlgebra;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_ref_frame
                     matrix_of_functions
                     find_euler_angles
                     flatten
                     matrix_product
                     mult_matrix_product
                     matrix_sum
                     matrix_sub
                     normalize
                     reshape
                     scalar_multipl
                     switch_ref_frame
                     translation
                     transpose
                     vector_cross
                     vector_length
                     x_axis_rotation
                     y_axis_rotation
                     z_axis_rotation );

use Carp qw( confess );
use Clone qw( clone );
use ForceField::Parameters;
use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------- Numeric linear algebra --------------------------- #

#
# Creates local reference frame for any three given atoms positions in cartesian
# coordinate system.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
# Output:
#     @local_ref_frame - Cartesian coordinates of points on x, y and z axis.
#

sub create_ref_frame
{
    my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) = @_;

    my @local_ref_frame;
    my $vector_length;

    # Let local z-axis be colinear to bond between mid and up atoms.
    $local_ref_frame[2][0] = $up_atom_coord->[0] - $mid_atom_coord->[0];
    $local_ref_frame[2][1] = $up_atom_coord->[1] - $mid_atom_coord->[1];
    $local_ref_frame[2][2] = $up_atom_coord->[2] - $mid_atom_coord->[2];

    # Let local x-axis be perpendicular to mid-up and mid-side bonds.
    $local_ref_frame[0][0] =
        ( $side_atom_coord->[1] - $mid_atom_coord->[1] ) * $local_ref_frame[2][2] -
        ( $side_atom_coord->[2] - $mid_atom_coord->[2] ) * $local_ref_frame[2][1];
    $local_ref_frame[0][1] =
      - ( $side_atom_coord->[0] - $mid_atom_coord->[0] ) * $local_ref_frame[2][2] +
        ( $side_atom_coord->[2] - $mid_atom_coord->[2] ) * $local_ref_frame[2][0];
    $local_ref_frame[0][2] =
        ( $side_atom_coord->[0] - $mid_atom_coord->[0] ) * $local_ref_frame[2][1] -
        ( $side_atom_coord->[1] - $mid_atom_coord->[1] ) * $local_ref_frame[2][0];

    # Let local y-axis be in the same plane as mid-up and mid-side bonds.
    $local_ref_frame[1][0] =
        $local_ref_frame[2][1] * $local_ref_frame[0][2] -
        $local_ref_frame[2][2] * $local_ref_frame[0][1];
    $local_ref_frame[1][1] =
      - $local_ref_frame[2][0] * $local_ref_frame[0][2] +
        $local_ref_frame[2][2] * $local_ref_frame[0][0];
    $local_ref_frame[1][2] =
        $local_ref_frame[2][0] * $local_ref_frame[0][1] -
        $local_ref_frame[2][1] * $local_ref_frame[0][0];

    # Normalizes all vectors to unit vectors.
    $vector_length =
        vector_length( [ [ $local_ref_frame[2][0] ],
                         [ $local_ref_frame[2][1] ],
                         [ $local_ref_frame[2][2] ] ] );
    $local_ref_frame[2][0] = $local_ref_frame[2][0] / $vector_length;
    $local_ref_frame[2][1] = $local_ref_frame[2][1] / $vector_length;
    $local_ref_frame[2][2] = $local_ref_frame[2][2] / $vector_length;
    $vector_length =
        vector_length( [ [ $local_ref_frame[0][0] ],
                         [ $local_ref_frame[0][1] ],
                         [ $local_ref_frame[0][2] ] ] );
    $local_ref_frame[0][0] = $local_ref_frame[0][0] / $vector_length;
    $local_ref_frame[0][1] = $local_ref_frame[0][1] / $vector_length;
    $local_ref_frame[0][2] = $local_ref_frame[0][2] / $vector_length;
    $vector_length =
        vector_length( [ [ $local_ref_frame[1][0] ],
                         [ $local_ref_frame[1][1] ],
                         [ $local_ref_frame[1][2] ] ] );
    $local_ref_frame[1][0] = $local_ref_frame[1][0] / $vector_length;
    $local_ref_frame[1][1] = $local_ref_frame[1][1] / $vector_length;
    $local_ref_frame[1][2] = $local_ref_frame[1][2] / $vector_length;

    return \@local_ref_frame;
}

#
# Function calculates Euler rotational angles (alpha, beta, gamma) that are used
# to transform global reference frame to chosen one.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
# Output:
#     euler angles (alpha, beta, gamma) in radians.
#

sub find_euler_angles
{
    my ( $parameters, $mid_atom_coord, $up_atom_coord, $side_atom_coord ) = @_;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $epsilon = $parameters->{'_[local]_constants'}{'epsilon'};

    my $alpha_rad;
    my $beta_rad;
    my $gamma_rad;

    my $z_axis_in_xy_plane;

    my $local_ref_frame =
        create_ref_frame( $mid_atom_coord, $up_atom_coord, $side_atom_coord );

    # Projects local z-axis to global xy-plane.
    $z_axis_in_xy_plane =
        sqrt( $local_ref_frame->[2][0] * $local_ref_frame->[2][0] +
              $local_ref_frame->[2][1] * $local_ref_frame->[2][1] );

    if( $z_axis_in_xy_plane > $epsilon ) {
        $alpha_rad =
            atan2  $local_ref_frame->[1][0] * $local_ref_frame->[2][1] -
                   $local_ref_frame->[1][1] * $local_ref_frame->[2][0],
                   $local_ref_frame->[0][0] * $local_ref_frame->[2][1] -
                   $local_ref_frame->[0][1] * $local_ref_frame->[2][0];
        $beta_rad = atan2 $z_axis_in_xy_plane, $local_ref_frame->[2][2];
        $gamma_rad =
            - atan2 - $local_ref_frame->[2][0], $local_ref_frame->[2][1];
    } else {
        $alpha_rad = 0.;
        $beta_rad = ( $local_ref_frame->[2][2] > 0. ) ? 0. : $pi;
        $gamma_rad = - atan2 $local_ref_frame->[0][1], $local_ref_frame->[0][0];
    }

    return [ $alpha_rad, $beta_rad, $gamma_rad ];
}

#
# Switches between global and local frames of reference.
# Input:
#     ${mid,up,side}_atom_coord - array of three atom coordinates in x, y, z
#     form, that frame of reference will be transformed to.
#     $switch_to_local - string that switches frame of reference to local (if
#     'local') and global (if 'global').
# Output:
#     $ref_frame_switch - 4x4 transformation matrix.
#

sub switch_ref_frame
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $switch_ref_to ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha_rad, $beta_rad, $gamma_rad ) =
        @{ find_euler_angles( $parameters,
                              $mid_atom_coord,
                              $up_atom_coord,
                              $side_atom_coord ) };

    # Depending on the option switch_to_local,
    my $ref_frame_switch;

    if( $switch_ref_to eq 'local' ) {
        $ref_frame_switch =
            mult_matrix_product( [ z_axis_rotation( $alpha_rad ),
                                   x_axis_rotation( $beta_rad ),
                                   z_axis_rotation( $gamma_rad ),
                                   translation( ( - $mid_atom_coord->[0],
                                                  - $mid_atom_coord->[1],
                                                  - $mid_atom_coord->[2] ) ) ] );
    } elsif( $switch_ref_to eq 'global' ) {
        $ref_frame_switch =
            mult_matrix_product( [ translation( ( $mid_atom_coord->[0],
                                                  $mid_atom_coord->[1],
                                                  $mid_atom_coord->[2] ) ),
                                   z_axis_rotation( - $gamma_rad ),
                                   x_axis_rotation( - $beta_rad ),
                                   z_axis_rotation( - $alpha_rad ) ] );
    } else {
        confess 'must choose $switch_to_global value between \'local\' and ' .
                '\'global\'';
    }

    return $ref_frame_switch;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around x-axis.
# Input:
#     $angle_rad - angle of rotation in radians.
# Output:
#     @rot_matrix_x - 4x4 matrix.
#

sub x_axis_rotation
{
    my ( $angle_rad ) = @_;

    my @rot_matrix_x =
        ( [ 1, 0, 0, 0 ],
          [ 0, cos( $angle_rad ), - sin( $angle_rad ), 0 ],
          [ 0, sin( $angle_rad ), cos( $angle_rad ), 0 ],
          [ 0, 0, 0, 1 ], );

    return \@rot_matrix_x;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around y-axis.
# Input:
#     $angle_rad - angle of rotation in radians.
# Output:
#     @rot_matrix_y - 4x4 matrix.
#

sub y_axis_rotation
{
    my ( $angle_rad ) = @_;

    my @rot_matrix_y =
        ( [ cos( $angle_rad ), 0, sin( $angle_rad ), 0 ],
          [ 0, 1, 0, 0 ],
          [ - sin( $angle_rad ), 0, cos( $angle_rad ), 0 ],
          [ 0, 0, 0, 1 ], );

    return \@rot_matrix_y;
}

#
# Creates 4x4 matrix that can rotate 4x4 matrices around z-axis.
# Input:
#     $angle_rad - angle of rotation in radians.
# Output:
#     @rot_matrix_z - 4x4 matrix.
#

sub z_axis_rotation
{
    my ( $angle_rad ) = @_;

    my @rot_matrix_z =
        ( [ cos( $angle_rad ), - sin( $angle_rad ), 0, 0 ],
          [ sin( $angle_rad ), cos( $angle_rad ), 0, 0 ],
          [ 0, 0, 1, 0 ],
          [ 0, 0, 0, 1 ], );

    return \@rot_matrix_z;
}

#
# Creates 4x4 matrix that can translate 4x4 matrices.
# Input:
#     @transl_coord - 3x1 matrix of x, y, z coordinates for corresponding
#     displacement.
# Output:
#     @transl_matrix - 4x4 matrix.
#

sub translation
{
    my @transl_coord = @_;

    my @transl_matrix =
        ( [ 1, 0, 0, $transl_coord[0] ],
          [ 0, 1, 0, $transl_coord[1] ],
          [ 0, 0, 1, $transl_coord[2] ],
          [ 0, 0, 0, 1 ], );

    return \@transl_matrix;
}

#
# Reshapes matrices according to given matrix dimensions.
# Input:
#     $element_list - array of all elements in matrix or matrices.
#     $dimensions - desired dimensions (m x n) of matrix or matrices in an array
#     form.
# Output:
#     @reshaped_matrices - array of matrix or matrices of desired dimensions.
#

sub reshape
{
    my ( $element_list, $dimensions ) = @_;

    # Checks, if there are enough elements for given list of dimensions.
    my $length_by_dimensions = 0;
    for( my $i = 0; $i < scalar @{ $dimensions }; $i += 2 ) {
        $length_by_dimensions += $dimensions->[$i] * $dimensions->[$i+1];
    }

    if( scalar( @{ $element_list } ) != $length_by_dimensions ) {
        confess 'there are not enough elements or dimensions';
    }

    # Generates matrices.
    my @matrices;
    for( my $i = 0; $i < scalar @{ $dimensions }; $i += 2 ) {
        my $rows = $dimensions->[$i];
        my $cols = $dimensions->[$i+1];

        push @matrices, [];

        foreach( 0..$rows-1 ) {
            push @{ $matrices[-1] }, [ splice @{ $element_list }, 0, $cols ];
        }
    }

    return \@matrices;
}

#
# Flattens all matrices to a list of elements.
# Input:
#     $matrices - array of matrix or matrices.
# Output:
#     @element_list - array of all elements in matrix or matrices.
#

sub flatten
{
    my ( $matrices ) = @_;

    my @element_list;

    for my $matrix ( @{ $matrices } ) {
    for my $row ( @{ $matrix } ) {
    for my $element ( @{ $row } ) {
        push @element_list, $element;
    } } };

    return \@element_list;
}

#
# Calculates vector length.
# Input:
#     $vector - 3x1 (if 4x1, last row is ignored) matrix.
# Output:
#     $vector_length - vector length.
#

sub vector_length
{
    my ( $vector ) = @_;

    my $vector_length = sqrt( $vector->[0][0] ** 2 +
                              $vector->[1][0] ** 2 +
                              $vector->[2][0] ** 2 );

    return $vector_length;
}

#
# Produces cross product of two 3D vectors.
# Input:
#     $left_matrix - left 3D vector where units are i, j, k.
#     $right_matrix - right 3D vector where units are i, j, k.
# Output:
#     @cross_product - cross product.
#

sub vector_cross
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @cross_product =
        ( $left_matrix->[1] * $right_matrix->[2] -
          $left_matrix->[2] * $right_matrix->[1],
        - $left_matrix->[0] * $right_matrix->[2] +
          $left_matrix->[2] * $right_matrix->[0],
          $left_matrix->[0] * $right_matrix->[1] -
          $left_matrix->[1] * $right_matrix->[0], );

    return \@cross_product;
}

#
# Converts vectors to unit vectors (normalizes).
# Input:
#     $vector - vector where units are i, j, k.
# Output:
#     @vector_normalized - normalized vector.
#

sub normalize
{
    my ( $vector ) = @_;

    my $vector_length = vector_length( $vector );

    my @vector_normalized;
    for my $row ( @{ $vector } ) {
        push @vector_normalized, [ $row->[0] / $vector_length ];
    }

    return \@vector_normalized;
}

#
# Adds two matrices.
# Input:
#     $left_matrix - left matrix.
#     $right_matrix - left matrix.
# Output:
#     $matrix_sum - sum of matrices (order of input variables is unimportant).
#

sub matrix_sum
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @matrix_sum;

    for my $i ( 0..$#{ $left_matrix } ) {
        for my $j ( 0..$#{ $left_matrix->[$i] } ) {
            if(!defined $right_matrix->[$i] || !defined $right_matrix->[$i][$j]){
                confess "there is no ($i, $j) element in the right matrix";
            }

            $matrix_sum[$i][$j] =
                $left_matrix->[$i][$j] + $right_matrix->[$i][$j];
        }
    }

    return \@matrix_sum;
}

#
# Subtracts two matrices.
# Input:
#     $left_matrix - left matrix.
#     $right_matrix - left matrix.
# Output:
#     $matrix_sub - subtraction of matrices (order of input variables is
#     unimportant).
#

sub matrix_sub
{
    my ( $left_matrix, $right_matrix ) = @_;

    my @matrix_sub;

    for my $i ( 0..$#{ $left_matrix } ) {
        for my $j ( 0..$#{ $left_matrix->[$i] } ) {
            if(!defined $right_matrix->[$i] || !defined $right_matrix->[$i][$j]){
                confess "there is no ($i, $j) element in the right matrix";
            }

            $matrix_sub[$i][$j] =
                $left_matrix->[$i][$j] - $right_matrix->[$i][$j];
        }
    }

    return \@matrix_sub;
}

#
# Multiplies matrix by scalar.
# Input:
#     $matrix - matrix.
#     $scalar - scalar value.
# Output:
#     @matrix - matrix multiplied by scalar.
#

sub scalar_multipl
{
    my ( $matrix, $scalar ) = @_;

    my @matrix_multipl;

    for my $i ( 0..$#{ $matrix } ) {
        for my $j ( 0..$#{ $matrix->[$i] } ) {
            $matrix_multipl[$i][$j] =
                $matrix->[$i][$j] * $scalar;
        }
    }

    return \@matrix_multipl;
}

# ------------------------- Matrices with functions --------------------------- #

sub matrix_of_functions
{
    my ( $function_ref, $rows, $cols ) = @_;
    my @matrix_of_functions = ();
    for my $row ( 0..$rows-1 ) {
        for my $col ( 0..$cols-1 ) {
            $matrix_of_functions[$row][$col] = $function_ref;
        }
    }
    return \@matrix_of_functions;
}

# ------------------------- Symbolic linear algebra --------------------------- #

#
# Performs basic linear algebra on symbolic expressions.
#
# Example of rotation along z-axis by chi angle:
#
#      / cos(chi) -sin(chi) 0 \   / x \   / x * cos(chi) + y * sin(chi) \
#      | sin(chi)  cos(chi) 0 | * | y | = | x * sin(chi) + y * cos(chi) |
#      \    0         0     1 /   \ z /   \              z              /
#

#
# Transposes matrix.
# Input:
#     $matrix - array representing matrix.
# Output:
#     @transposed_matrix - transposed matrix.
#

sub transpose
{
    my ( $matrix ) = @_;
    my @matrix = @{ $matrix };

    my @transposed_matrix;

    for my $row ( 0..$#matrix ) {
        for my $col ( 0..$#{ $matrix[$row] } ) {
            $transposed_matrix[$col][$row] = $matrix[$row][$col];
        }
    }

    return \@transposed_matrix;
}

#
# Calculates matrix product of two matrices.
# Input:
#     ${left, right}_matrix - matrices.
#     $symbol_values - values of the unknown variable(-s).
# Output:
#     @matrix_product - matrix product.
#

sub matrix_product
{
    my ( $left_matrix, $right_matrix, $symbol_values ) = @_;

    my %symbol_values;
    if( defined $symbol_values ) { %symbol_values = %{ $symbol_values } };

    # Checks for analytical functions and evaluates if the values are present.
    if( ref $left_matrix eq 'Symbolic' ) {
        if( $left_matrix->{'is_evaluated'} ) {
            $left_matrix = $left_matrix->{'matrix'};
        } else {
            my @symbol_values =
                map { $symbol_values->{$_} } @{ $left_matrix->{'symbols'} };
            if( defined $symbol_values[-1] &&
                $#symbol_values eq $#{ $left_matrix->{'symbols'} } ) {
                $left_matrix->evaluate( $symbol_values );
                $left_matrix = $left_matrix->{'matrix'};
            }
        }
    }

    if( ref $right_matrix eq 'Symbolic' ) {
        if( $right_matrix->{'is_evaluated'} ) {
            $right_matrix = $right_matrix->{'matrix'};
        } else {
            my @symbol_values =
                map { $symbol_values->{$_} } @{ $right_matrix->{'symbols'} };
            if( defined $symbol_values[-1] &&
                $#symbol_values eq $#{ $right_matrix->{'symbols'} } ) {
                $right_matrix->evaluate( $symbol_values );
                $right_matrix = $right_matrix->{'matrix'};
            }
        }
    }

    # Returns both matrices if any of two matrices are not evaluated.
    if( ref $left_matrix eq 'Symbolic' || ref $right_matrix eq 'Symbolic' ) {
        return $left_matrix, $right_matrix;
    }

    # Notifies error, when the column number of left matrix does not equal the
    # row number of the right matrix.
    if( scalar @{ transpose( $left_matrix ) } != scalar @{ $right_matrix } ) {
        confess { type => 'DimensionError',
              message => 'A row number of a left matrix is NOT equal ' .
                         "to the column\nnumber of the right matrix.", };
    }

    # Calculates matrix product of two matrices.
    my @matrix_product;
    for(my $left_row = 0; $left_row <= $#{ $left_matrix }; $left_row++){
    for(my $right_col = 0; $right_col <= $#{ $right_matrix->[0] }; $right_col++){
    for(my $right_row = 0; $right_row <= $#{ $right_matrix }; $right_row++){
        if( ! defined $matrix_product[$left_row][$right_col] ) {
            $matrix_product[$left_row][$right_col] =
                $left_matrix->[$left_row][$right_row] *
                $right_matrix->[$right_row][$right_col];
        } else {
            $matrix_product[$left_row][$right_col] +=
                $left_matrix->[$left_row][$right_row] *
                $right_matrix->[$right_row][$right_col];
        }
    } } }

    return \@matrix_product;
}

#
# Calculates matrix product of list of any size of matrices.
# Input:
#     $matrices - list of matrices.
#     $symbol_values - values of the unknown variable(-s).
# Output:
#     @mult_matrix_product - matrix product.
#

sub mult_matrix_product
{
    my ( $matrices, $symbol_values ) = @_;

    my @matrices = @{ clone( $matrices ) };

    # Multiplies matrices from left to right.
    my @mult_matrix_product;

    # Only evaluates variables if there is only one matrix in @matrices.
    if( scalar @matrices == 1 ) {
        my $matrix = $matrices[0];
        # Checks for analytical functions and evaluates if the values are
        # present.
        if( ref $matrix eq 'Symbolic' ) {
            if( $matrix->{'is_evaluated'} ) {
                $matrix = $matrix->{'matrix'};
            } else {
                my @symbol_values =
                    map { $symbol_values->{$_} } @{ $matrix->{'symbols'} };
                if( defined $symbol_values[-1] &&
                    $#symbol_values eq $#{ $matrix->{'symbols'} } ) {
                    $matrix->evaluate( $symbol_values );
                    $matrix = $matrix->{'matrix'};
                }
            }
        }

        push @mult_matrix_product, $matrix;

    } else {
        for( my $i = $#matrices; $i >= 1; $i-- ) {
            if( $i == $#matrices ) {
                eval {
                    push @mult_matrix_product,
                        matrix_product( $matrices[$i-1],
                                        $matrices[$i],
                                        $symbol_values );
                } or do {
                    confess ${ @ }->{message};
                };
            } else {
                eval {
                    unshift @mult_matrix_product,
                            matrix_product( $matrices[$i-1],
                                            $mult_matrix_product[0],
                                            $symbol_values );
                    splice @mult_matrix_product, 1, 1;
                } or do {
                    confess ${ @ }->{message};
                };
            }
        }
    }

    return \@mult_matrix_product;
}

1;
