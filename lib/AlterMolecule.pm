package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending
                     bond_stretching
                     bond_torsion );

use LinearAlgebra qw( mult_matrix_product
                      switch_ref_frame );

# ----------------------- Molecule alteration matrices ------------------------ #

#
# Makes a rotational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_name - name of the dihedral angle.
# Output:
#     $rot_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined dihedral angle variables).
#

sub bond_torsion
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $angle_name ) = @_;

    # Rotation matrix around the bond.
    my @rot_matrix =
	( [ "cos(\$${angle_name})", "-sin(\$${angle_name})", 0, 0 ],
	  [ "sin(\$${angle_name})",  "cos(\$${angle_name})", 0, 0 ],
	  [ 0, 0, 1, 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
    	mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
    						    $up_atom_coord,
    						    $side_atom_coord,
    						    "global" ) },
    			       \@rot_matrix,
    			       @{ switch_ref_frame( $mid_atom_coord,
    						    $up_atom_coord,
    						    $side_atom_coord,
    						    "local" ) } ] );

    return $rot_matrix;
}

#
# Makes a translational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $length_name - name of the bond length variable.
# Output:
#     $transl_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined bond length variable).
#

sub bond_stretching
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $length_name ) = @_;

    # Translation of the coordinates of the bond.
    my @transl_matrix = ( [ 1, 0, 0, 0 ],
			  [ 0, 1, 0, 0 ],
			  [ 0, 0, 1, "\$${length_name}" ],
			  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $transl_matrix =
    	mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
						    $up_atom_coord,
						    $side_atom_coord,
						    "global" ) },
			       \@transl_matrix,
			       @{ switch_ref_frame( $mid_atom_coord,
						    $up_atom_coord,
						    $side_atom_coord,
						    "local" ) } ] );

    return $transl_matrix;
}

#
# Makes a transformation matrix for changing angles between two bonds.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_name_x - name of the bond angle that will be rotated in x-axis.
#     $angle_name_y - name of the bond angle that will be rotated in y-axis.
# Output:
#     $rot_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined angle variables).
#

sub angle_bending
{
    my ( $mid_atom_coord,
	 $up_atom_coord,
	 $side_atom_coord,
	 $angle_name_x,
	 $angle_name_y, ) = @_;

    # Bond angle matrices that rotates along x and y axes.
    my @rot_matrix_x =
	( [ 1, 0, 0, 0 ],
	  [ 0, "cos(\$${angle_name_x})", "-sin(\$${angle_name_x})", 0 ],
	  [ 0, "sin(\$${angle_name_x})", "cos(\$${angle_name_x})", 0 ],
	  [ 0, 0, 0, 1 ] );
    my @rot_matrix_y =
	( [ "cos(\$${angle_name_y})", 0, "sin(\$${angle_name_y})", 0 ],
	  [ 0, 1, 0, 0 ],
	  [ "-sin(\$${angle_name_y})", 0, "cos(\$${angle_name_y})", 0 ],
	  [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
	mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
						    $up_atom_coord,
						    $side_atom_coord,
						    "global" ) },
			       \@rot_matrix_y,
			       \@rot_matrix_x,
			       @{ switch_ref_frame( $mid_atom_coord,
						    $up_atom_coord,
						    $side_atom_coord,
						    "local" ) } ] );

    return $rot_matrix;
}

1;
