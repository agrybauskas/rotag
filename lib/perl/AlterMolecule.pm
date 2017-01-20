package AlterMolecule;

use strict;
use warnings;

use lib qw( ./ );
use LinearAlgebra;

use Data::Dumper;

#
# Makes a rotational transformation matrix of the bond.
# Input  ():
# Output ():
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
	[ [ 1, 0, 0, ( -1 ) * $mid_atom_coord->[0] ],
	  [ 0, 1, 0, ( -1 ) * $mid_atom_coord->[1] ],
	  [ 0, 0, 1, ( -1 ) * $mid_atom_coord->[2] ],
	  [ 0, 0, 0,                   1           ] ];

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
			 [ 0, cos( $beta ), ( -1 ) * sin( $beta ), 0 ],
			 [ 0, sin( $beta ), cos( $beta ), 0 ],
			 [ 0, 0, 0, 1 ] );

    my @rot_matrix_y = ( [ cos( $gamma ), 0, sin( $gamma ), 0 ],
			 [ 0, 1, 0, 0 ],
			 [ ( -1 ) * sin( $gamma ), 0, cos( $gamma ), 0 ],
			 [ 0, 0, 0, 1 ]);

    my @rot_matrix_z = ( [ cos( $alpha ), (-1) * sin( $alpha ), 0, 0 ],
			 [ sin( $alpha ), cos( $alpha ), 0, 0 ],
			 [ 0, 0, 1, 0 ],
			 [ 0, 0, 0, 1 ] );

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
