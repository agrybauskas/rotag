package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending
                     bond_stretching
                     bond_torsion );

use lib qw( ./ );
use LinearAlgebra qw( find_euler_angles
                      matrix_product
                      switch_ref_frame );

# ------------------ Molecule structure alteration algorithms ----------------- #

#
# Changes structure of the molecule by rotating along dihedral angles, changing
# length and angle of the bonds.
#

#
# Makes a rotational transformation matrix of the bond.
# Input  (4 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side and symbol describing dihedral angle.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined dihedral angle).
#

sub bond_torsion
{
    my ( $angle_symbol,
	 $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix around the bond.
    my @bond_rot_matrix =
	( [ "cos(${angle_symbol})", "-sin(${angle_symbol})", 0, 0 ],
	  [ "sin(${angle_symbol})",  "cos(${angle_symbol})", 0, 0 ],
	  [ 0, 0, 1, 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
    	&matrix_product( &switch_ref_frame( "global",
					    $mid_atom_coord,
					    $up_atom_coord,
					    $side_atom_coord ),
			 \@bond_rot_matrix,
			 &switch_ref_frame( "local",
					    $mid_atom_coord,
					    $up_atom_coord,
					    $side_atom_coord ) );

    return $rot_matrix;
}

#
# Changes the length of the bond.
# Input  (4 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side and symbol for bond length.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined bond length variable: bond_length).
#

sub bond_stretching
{
    my ( $length_symbol,
	 $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	find_euler_angles( @$mid_atom_coord,
			   @$up_atom_coord,
			   @$side_atom_coord );

    # Translation of the bond.
    my @bond_len_matrix = ( [ 1, 0, 0, 0 ],
			    [ 0, 1, 0, 0 ],
			    [ 0, 0, 1, "${length_symbol}" ],
			    [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $transl_matrix =
	matrix_product( &switch_ref_frame( "global",
					   $mid_atom_coord,
					   $up_atom_coord,
					   $side_atom_coord ),
			\@bond_len_matrix,
			&switch_ref_frame( "local",
					   $mid_atom_coord,
					   $up_atom_coord,
					   $side_atom_coord ) );

    return $transl_matrix;
}

#
# Changes the length of the bond.
# Input  (3 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined bond angle variables: theta and psi).

sub angle_bending
{
    my ( $first_angle_symbol,
	 $second_angle_symbol,
	 $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	find_euler_angles( @$mid_atom_coord,
			   @$up_atom_coord,
			   @$side_atom_coord );

    # Bond angle matrices.
    my @rot_x_matrix =
	( [ 1, 0, 0, 0 ],
	  [ 0, "cos(${first_angle_symbol})", "-sin(${first_angle_symbol})", 0 ],
	  [ 0, "sin(${first_angle_symbol})", "cos(${first_angle_symbol})", 0 ],
	  [ 0, 0, 0, 1 ] );
    my @rot_y_matrix =
	( [ "cos(${second_angle_symbol})", 0, "sin(${second_angle_symbol})", 0 ],
	  [ 0, 1, 0, 0 ],
	  [ "-sin(${second_angle_symbol})", 0, "cos(${second_angle_symbol})", 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $angle_matrix =
	matrix_product( &switch_ref_frame( "global",
					   $mid_atom_coord,
					   $up_atom_coord,
					   $side_atom_coord ),
			\@rot_y_matrix,
			\@rot_x_matrix,
			&switch_ref_frame( "local",
					   $mid_atom_coord,
					   $up_atom_coord,
					   $side_atom_coord ) );

    return $angle_matrix;
}

1;
