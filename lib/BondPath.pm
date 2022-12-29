package BondPath;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my $self = { 'value' => undef };
    return bless $self, $class;
}

# ------------------------- Struct-related functions -------------------------- #

sub stretchable_bonds
{

}

sub bendable_angles
{

}

sub rotatable_bonds
{

}

1;
