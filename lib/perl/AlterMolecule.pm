package AlterMolecule;

use strict;
use warnings;

use lib qw( ./ );
use LinearAlgebra;

use Data::Dumper;

# ------------------ Molecule structure alteration algorithms ----------------- #

#
# Block of code that has functions changing structure of molecules by rotating
# along dihedral angles, changing length and angle of bonds.
#

#
# Makes a rotational transformation matrix of the bond.
# Input  (4 arg): cartesian coordinates in array form that defines user-selected
#                 mid, up, side and target atoms' coordinates.
# Output (1 arg): matrix defining coordinates in symbolic mathematical form
#                 (with undefined dihedral angle chi).
#

sub rotate_bond
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $target_atom_coord ) = @_;

    # Transformations that transforms global reference frame to local.
    # Negative translation matrix. Places mid-atom to 
    my @neg_transl_matrix =
	( [ 1, 0, 0, - $mid_atom_coord->[0] ],
	  [ 0, 1, 0, - $mid_atom_coord->[1] ],
	  [ 0, 0, 1, - $mid_atom_coord->[2] ],
	  [ 0, 0, 0,            1           ] );

    # Positive translation matrix. Brings back to original mid-atom position.
    my @pos_transl_matrix =
	( [ 1, 0, 0, $mid_atom_coord->[0] ],
	  [ 0, 1, 0, $mid_atom_coord->[1] ],
	  [ 0, 0, 1, $mid_atom_coord->[2] ],
	  [ 0, 0, 0, 1 ] );

    # Rotation matrix to coordinating global reference frame properly.
    # Finding Euler angles necessary for rotation matrix.
    my ( $alpha, $beta, $gamma ) =
	LinearAlgebra::find_euler_angles( @$mid_atom_coord,
					  @$up_atom_coord,
					  @$side_atom_coord );

    my @rot_matrix_x = ( [ 1, 0, 0, 0 ],
			 [ 0, cos( $beta ), - sin( $beta ), 0 ],
			 [ 0, sin( $beta ),   cos( $beta ), 0 ],
			 [ 0, 0, 0, 1 ] );

    my @rot_matrix_y = ( [   cos( $gamma ), 0, sin( $gamma ), 0 ],
			 [ 0, 1, 0, 0 ],
			 [ - sin( $gamma ), 0, cos( $gamma ), 0 ],
			 [ 0, 0, 0, 1 ] );

    my @rot_matrix_z = ( [ cos( $alpha ), - sin( $alpha ), 0, 0 ],
			 [ sin( $alpha ),   cos( $alpha ), 0, 0 ],
			 [ 0, 0, 1, 0 ],
			 [ 0, 0, 0, 1 ] );

    # Rotation matrix around the bond.
    my @rot_matrix_bond = ( [ 'cos($chi)', '-sin($chi)', 0, 0 ],
			    [ 'sin($chi)',  'cos($chi)', 0, 0 ],
			    [ 0, 0, 1, 0 ],
			    [ 0, 0, 0, 1 ] );

    # Converting target atom coordinates from 3x1 to 4x1 form so,
    # it could be multiplied by 4x4 matrix.
    my @target_atom_coord;

    foreach( @$target_atom_coord ) {
    	push( @target_atom_coord, [ $_ ] );
    }

    push( @target_atom_coord, [ 1 ] );

    # Incorporating symbol for dihedral chi angle.
    my @symbols = [ "chi", "i" ];

    # Multiplying multiple matrices to get a final form.
    my @rot_matrix;

    @rot_matrix = LinearAlgebra::mult_matrix_product(
    	@symbols,
    	\@pos_transl_matrix,
    	LinearAlgebra::transpose( \@rot_matrix_y ),
    	LinearAlgebra::transpose( \@rot_matrix_x ),
    	LinearAlgebra::transpose( \@rot_matrix_z ),
	\@rot_matrix_bond,
    	\@rot_matrix_z,
    	\@rot_matrix_x,
    	\@rot_matrix_y,
	\@neg_transl_matrix,
    	\@target_atom_coord );

    return \@rot_matrix;
}

#
# Changes the length of the bond.
# Input  ():
# Ouptut ():
#

sub change_bond_len
{

}

#
# Changes the length of the bond.
# Input  ():
# Ouptut ():
#

sub change_bond_angle
{

}

1;
