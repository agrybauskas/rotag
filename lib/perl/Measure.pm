package Measure;

use strict;
use warnings;

use Math::Trig;

use lib "./";
use LinearAlgebra;

use Data::Dumper;

# ----------------------------- Molecule parameters --------------------------- #

#
# Calculates various parameters that describe molecule or atoms, such as, bond
# length, dihedral angle, torsion angle, RMSD and etc.
#

#
# Calculates bond length of given two atoms.
# Input  (2 arg): matrices of x,y,z coordinates of two atoms.
# Output (1 arg): length of the bond in angstroms.
#

sub bond_length
{
    my @atom_coord = @_;

    my $bond_length =
	sqrt( ( $atom_coord[1][0] - $atom_coord[0][0] )**2
	    + ( $atom_coord[1][1] - $atom_coord[0][1] )**2
	    + ( $atom_coord[1][2] - $atom_coord[0][2] )**2 );

    return $bond_length;
}

#
# Calculates angle between three atoms.
# Input  (3 arg): matrices of x,y,z coordinates of three atoms.
# Output (1 arg): angle in radians.
#

sub bond_angle
{
    my $atom_coord = shift;
    my @atom_coord = @$atom_coord;

    my $bond_angle;

    # Angle between three atoms (in radians) in 3-D space can be calculated by
    # the formula:
    #                            ->   ->      ->         ->
    #            theta = arccos( AB * BC / || AB || * || BC || )

    # This formula is applied to atoms where vectors are the substraction of
    # coordinates of two atoms. Suppose, one of the side atom is A, B - middle
    # and C - remaining atom. Order of side atoms is irrelevant.
    my @vector_ab = ( $atom_coord[0][0] - $atom_coord[1][0],
		      $atom_coord[0][1] - $atom_coord[1][1],
		      $atom_coord[0][2] - $atom_coord[1][2] );
    my @vector_bc = ( $atom_coord[2][0] - $atom_coord[1][0],
		      $atom_coord[2][1] - $atom_coord[1][1],
		      $atom_coord[2][2] - $atom_coord[1][2] );

    my $vector_product = $vector_ab[0] * $vector_bc[0] +
	               + $vector_ab[1] * $vector_bc[1]
	               + $vector_ab[2] * $vector_bc[2];

    my $length_ab = sqrt( $vector_ab[0]**2
			+ $vector_ab[1]**2
			+ $vector_ab[2]**2 );
    my $length_bc = sqrt( $vector_bc[0]**2
			+ $vector_bc[1]**2
			+ $vector_bc[2]**2 );

    $bond_angle = acos( $vector_product / ( $length_ab * $length_bc ) );

    return $bond_angle;
}

#
# Calculates dihedral angle of four given atoms.
# Input  (4 arg): matrices of x,y,z coordinates of four atoms.
# Output (1 arg): dihedral angle in radians.
#

sub dihedral_angle
{
    my $atom_coord = shift;
    my @atom_coord = @$atom_coord;

    my $dihedral_angle;

    # Dihedral angle is calculated by, first, calculating vectors that
    # lie on the surface described by three points using equation of the
    # plane:
    #
    #                          ax * by * cy = d
    #
    # where x, y and z are values for cartesian coordinates and a, b, c
    # d - constants describing plane orientation in space.
    #
    # Dihedral angle is calculated by using constants a, b and c, that
    # describes two planes (plane_left, plane_right):
    #
    # psi = arccos( a_1 * a_2 + b_1 * b_2 + c_1 * c_2
    #             / sqrt( a_1^2 + b_1^2 + c_1^2 )
    #             * sqrt( a_2^2 + b_2^2 + c_2^2 ) )
    #

    # Calculates coefficients (a, b, c) for both normal vectors to the surfaces
    # described by three atom-coordinates each and sharing two.
    my @vector_left = map { $atom_coord[0][$_] - $atom_coord[1][$_] }
                          ( 0..$#{ $atom_coord[0] } );

    my @vector_right = map { $atom_coord[2][$_] - $atom_coord[1][$_] }
                           ( 0..$#{ $atom_coord[0] } );

    my ( $a_1, $b_1, $c_1) =
	( $vector_left[1] * $vector_right[2]
	- $vector_left[2] * $vector_right[1],
        - $vector_left[0] * $vector_right[2]
	+ $vector_left[2] * $vector_right[0],
	  $vector_left[1] * $vector_right[2]
	- $vector_left[2] * $vector_right[1] );

    @vector_left = map { $atom_coord[1][$_] - $atom_coord[2][$_] }
                       ( 0..$#{ $atom_coord[2] } );

    @vector_right = map { $atom_coord[3][$_] - $atom_coord[2][$_] }
                        ( 0..$#{ $atom_coord[2] } );

    my ( $a_2, $b_2, $c_2) =
	( $vector_left[1] * $vector_right[2]
	- $vector_left[2] * $vector_right[1],
        - $vector_left[0] * $vector_right[2]
	+ $vector_left[2] * $vector_right[0],
	  $vector_left[1] * $vector_right[2]
	- $vector_left[2] * $vector_right[1] );

    # Calculates dihedral angle in radians.
    $dihedral_angle = acos( ( $a_1 * $a_2 + $b_1 * $b_2 + $c_1 * $c_2 )
			    / ( sqrt( $a_1**2 + $b_1**2 + $c_1**2 )
			      * sqrt( $a_2**2 + $b_2**2 + $c_2**2 ) ) );

    return $dihedral_angle;
}

#
#
# Input  (  arg):
# Output (  arg):
#

sub rmsd
{

}

1;
