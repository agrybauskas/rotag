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

    #                  -> ->    ->
    # Creates vectors: a, b and c, that are translated to global reference frame.
    # Picture of vectors:
    #                                   ->  O ->
    #                                   b  /  c
    #                              -> CA---C
    #                              a /
    #                               N
    my $vector_a =
    	LinearAlgebra::matrix_sub( [ $atom_coord[1] ], [ $atom_coord[0] ] );
    my $vector_b =
    	LinearAlgebra::matrix_sub( [ $atom_coord[2] ], [ $atom_coord[1] ] );
    my $vector_c =
    	LinearAlgebra::matrix_sub( [ $atom_coord[3] ], [ $atom_coord[2] ] );

    #                                               ->    -> ->    ->
    # Creates normal vectors from the vector pairs: a and b, b and c.
    my $vector_cross_ab =
    	LinearAlgebra::vector_cross( @$vector_a, @$vector_b );
    my $vector_cross_bc =
    	LinearAlgebra::vector_cross( @$vector_b, @$vector_c );

    # Calculates length for each cross product.
    my $vector_length_ab = sqrt( $vector_cross_ab->[0]**2
			       + $vector_cross_ab->[1]**2
			       + $vector_cross_ab->[2]**2 );
    my $vector_length_bc = sqrt( $vector_cross_bc->[0]**2
			       + $vector_cross_bc->[1]**2
			       + $vector_cross_bc->[2]**2 );

    # Calculates normal vectors for each cross product of two vectors.
    my @normal_vector_ab = map { $_ / $vector_length_ab } @$vector_cross_ab;
    my @normal_vector_bc = map { $_ / $vector_length_bc } @$vector_cross_bc;
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
