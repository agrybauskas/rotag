package BondCombinations;

use strict;
use warnings;

use Combinatorics qw( permutation );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;

    my $self = {
        'order' => $args->{'order'},
        'collection' => $args->{'collection'},
        'cache' => undef,
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

sub get_names
{
    my ( $self ) = @_;
    return $self->{'order'};
}

sub get_values
{
    my ( $self, $names ) = @_;
    $names //= $self->{'order'};
    my $permuted_values = permutation(
        scalar( @{ $names } ),
        [],
        [ map { $self->{'collection'}{$_} } @{ $names } ],
        []
    );
    return $permuted_values;
}

1;
