package AlterMolecule;

use strict;
use warnings;

use lib qw( ./ );
use LinearAlgebra;

# ------------------ Molecule structure alteration algorithms ----------------- #

#
# Changes structure of the molecule by rotating along dihedral angles, changing
# length and angle of the bonds.
#

#
# Makes a rotational transformation matrix of the bond.
# Input  (3 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined dihedral angle chi).
#

sub bond_torsion
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix around the bond.
    my @bond_rot_matrix = ( [ 'cos($chi)', '-sin($chi)', 0, 0 ],
			    [ 'sin($chi)',  'cos($chi)', 0, 0 ],
			    [ 0, 0, 1, 0 ],
			    [ 0, 0, 0, 1 ] );

    # Incorporating symbol for dihedral chi angle.
    my @symbols = ( "chi", "i" );

    # Multiplying multiple matrices to get a final form.
    my @rot_matrix =
    	LinearAlgebra::mult_matrix_product(
    	    \@symbols,
    	    LinearAlgebra::switch_ref_frame( "global",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ),
    	    \@bond_rot_matrix,
    	    LinearAlgebra::switch_ref_frame( "local",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ) );

    return \@rot_matrix;
}

#
# Changes the length of the bond.
# Input  (3 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined bond length variable: bond_length).
#

sub bond_streching
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	LinearAlgebra::find_euler_angles( @$mid_atom_coord,
					  @$up_atom_coord,
					  @$side_atom_coord );

    # Translation of the bond.
    my @bond_len_matrix = ( [ 1, 0, 0, 0 ],
			    [ 0, 1, 0, 0 ],
			    [ 0, 0, 1, '$r' ],
			    [ 0, 0, 0, 1 ] );

    # Incorporating symbol for bond length - len.
    my @symbols = ( "r", "i" );

    # Multiplying multiple matrices to get a final form.
    my @transl_matrix =
	LinearAlgebra::mult_matrix_product(
	    \@symbols,
	    LinearAlgebra::switch_ref_frame( "global",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ),
	    \@bond_len_matrix,
	    LinearAlgebra::switch_ref_frame( "local",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ) );

    return \@transl_matrix;
}

#
# Changes the length of the bond.
# Input  (3 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined bond angle variables: theta and psi).

sub angle_bending
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord ) = @_;

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	LinearAlgebra::find_euler_angles( @$mid_atom_coord,
					  @$up_atom_coord,
					  @$side_atom_coord );

    # Bond angle matrices.
    my @rot_x_matrix = ( [ 1, 0, 0, 0 ],
			 [ 0, 'cos($theta)', '-sin($theta)', 0 ],
			 [ 0, 'sin($theta)', 'cos($theta)', 0 ],
			 [ 0, 0, 0, 1 ] );
    my @rot_y_matrix = ( [ 'cos($psi)', 0, 'sin($psi)', 0 ],
			 [ 0, 1, 0, 0 ],
			 [ '-sin($psi)', 0, 'cos($psi)', 0 ],
			 [ 0, 0, 0, 1 ] );

    # Incorporating symbol for bond length - len.
    my @symbols = ( "theta", "psi", "i" );

    # Multiplying multiple matrices to get a final form.
    my @angle_matrix =
	LinearAlgebra::mult_matrix_product(
	    \@symbols,
	    LinearAlgebra::switch_ref_frame( "global",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ),
	    \@rot_y_matrix,
	    \@rot_x_matrix,
	    LinearAlgebra::switch_ref_frame( "local",
    	    				     $mid_atom_coord,
    	    				     $up_atom_coord,
    	    				     $side_atom_coord ) );

    return \@angle_matrix;
}

1;
