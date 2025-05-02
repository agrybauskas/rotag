package BondCombinations;

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
        my ( $name ) = keys %{ $bond_parameter };
        $self->{'collection'}{$name} = $bond_parameter->{$name};
        push @{ $self->{'order'} }, $name;
    }
    return;
}

sub remove
{

}

# --------------------------------- Methods ----------------------------------- #

sub get_all_names
{
    my ( $self ) = @_;
    return $self->{'order'};
}

sub get_all_values
{

}

1;
