package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending bond_stretching bond_torsion );

use lib qw( ./ );
use LinearAlgebra qw( matrix_product switch_ref_frame );

# ------------------ Molecule structure alteration algorithms ----------------- #

#
# Changes structure of the molecule by rotating along dihedral angles, changing
# length and angle of the bonds.
#

#
# Makes a rotational transformation matrix of the bond.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_symbol - symbol describing dihedral angle.
# Output:
#     $rot_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined dihedral angle variable).
#

sub bond_torsion
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $angle_symbol ) = @_;

    # Rotation matrix around the bond.
    my @rot_matrix =
	( [ "cos(${angle_symbol})", "-sin(${angle_symbol})", 0, 0 ],
	  [ "sin(${angle_symbol})",  "cos(${angle_symbol})", 0, 0 ],
	  [ 0, 0, 1, 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
    	matrix_product( switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "global" ),
			\@rot_matrix,
			switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "local" ) );

    return $rot_matrix;
}

#
# Makes a translational transformation matrix of the bond.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_symbol - symbol describing bond length.
# Output:
#     $transl_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined bond length variable).
#

sub bond_stretching
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $length_symbol ) = @_;

    # Translation of the coordinates of the bond.
    my @transl_matrix = ( [ 1, 0, 0, 0 ],
			  [ 0, 1, 0, 0 ],
			  [ 0, 0, 1, "${length_symbol}" ],
			  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $transl_matrix =
    	matrix_product( switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "global" ),
    			\@transl_matrix,
    			switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "local" ) );

    return $transl_matrix;
}

#
# Makes a transformation matrix for changing angles between two bonds.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_symbol_x - symbol describing bond angle that will be rotated in
#     x-axis.
#     $angle_symbol_y - symbol describing bond angle that will be rotated in
#     y-axis.
# Output:
#     $angle_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined angle variables).
#

sub angle_bending
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $angle_symbol_x,
	 $angle_symbol_y, ) = @_;

    # Bond angle matrices that rotates along x and y axes.
    my @rot_matrix_x =
	( [ 1, 0, 0, 0 ],
	  [ 0, "cos(${angle_symbol_x})", "-sin(${angle_symbol_x})", 0 ],
	  [ 0, "sin(${angle_symbol_x})", "cos(${angle_symbol_x})", 0 ],
	  [ 0, 0, 0, 1 ] );
    my @rot_matrix_y =
	( [ "cos(${angle_symbol_y})", 0, "sin(${angle_symbol_y})", 0 ],
	  [ 0, 1, 0, 0 ],
	  [ "-sin(${angle_symbol_y})", 0, "cos(${angle_symbol_y})", 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $angle_matrix =
	matrix_product( switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "global" ),
			\@rot_matrix_y,
			\@rot_matrix_x,
			switch_ref_frame( $mid_atom_coord,
					  $up_atom_coord,
					  $side_atom_coord,
					  "local" ) );

    return $angle_matrix;
}

1;
