package BondCombIter;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;

    my $self = {
        'order' => undef,
        'collection' => $args->{'collection'},
    };

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub add
{
    my ( $self, $bond_parameters ) = @_;
}

sub remove
{

}

# --------------------------------- Methods ----------------------------------- #

sub next
{

}

sub prev
{

}

sub reset
{

}

sub is_end
{

}

sub get_all_names
{

}

sub get_all_values
{

}

1;
