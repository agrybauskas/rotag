package Measure;

use strict;
use warnings;

# ----------------------------- Molecule parameters --------------------------- #

#
# Calculates various parameters that describe molecule or atoms, such as, bond
# length, dihedral angle, torsion angle, RMSD and etc.
#

#
# 
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
# 
# Input  (3  arg): 
# Output (  arg): 
#

sub dihedral_angle
{
    
}

#
# 
# Input  (  arg): 
# Output (  arg): 
#

sub bond_angle
{
    
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
