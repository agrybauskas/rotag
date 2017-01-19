package AlterMolecule;

use strict;
use warnings;

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
    my $neg_transl_matrix =
	[ [ 1, 0, 0, - ( $mid_atom_coord->[0] ) ],
	  [ 0, 1, 0, - ( $mid_atom_coord->[1] ) ],
	  [ 0, 0, 1, - ( $mid_atom_coord->[2] ) ],
	  [ 0, 0, 0,              1             ] ];

    # Positive translation matrix. Brings back to original mid-atom position.
    my $pos_transl_matrix =
	[ [ 1, 0, 0, $mid_atom_coord->[0] ],
	  [ 0, 1, 0, $mid_atom_coord->[1] ],
	  [ 0, 0, 1, $mid_atom_coord->[2] ],
	  [ 0, 0, 0,          1           ] ];
    
    # my ( $mid_atom_x,  $mid_atom_y,  $mid_atom_z,
    #      $up_atom_x,   $up_atom_y,   $up_atom_z,
    #      $side_atom_x, $side_atom_y, $side_atom_z ) = @_    
}

#
# Changes the length of the bond.
# Input  ():
# Ouptut ():
#

sub change_bond_len
{

}
1;
