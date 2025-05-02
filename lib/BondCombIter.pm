package BondCombIter;

use strict;
use warnings;

use Clone qw( clone );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;

    my $self = {
        'order' => $args->{'order'},
        'collection' => $args->{'collection'},
    };

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub add
{
    my ( $self, $bond_parameters ) = @_;
    for my $bond_parameter ( @{ $bond_parameters } ) {
        # $self->{'collection'}{$bond_parameter} = clone $bond_parameters->{$_};
    }
    return;
}

sub remove
{

}

# --------------------------------- Methods ----------------------------------- #

sub get_all_names
{

}

sub get_all_values
{

}

1;
