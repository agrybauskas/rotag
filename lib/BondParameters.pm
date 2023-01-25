package BondParameters;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my $self = { 'parameters' => $args->{'parameters'} };
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_rotatable_bonds
{

}

sub set_stretchable_bonds
{

}

sub set_bendable_angles
{

}

# --------------------------------- Methods ----------------------------------- #

sub calculate_dihedral_angles
{

}

sub calculate_bond_lengths
{

}

sub calculate_bond_angles
{

}

1;
