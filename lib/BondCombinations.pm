package BondCombinations;

use strict;
use warnings;

use List::Util qw( any );

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
        if( ! any { $_ eq $name } @{ $self->{'order'} }  ) {
            push @{ $self->{'order'} }, $name;
        }
        for my $value ( @{ $bond_parameter->{$name} } ) {
            if( ! any { $_ eq $value } @{ $self->{'collection'}{$name} }  ) {
                push @{ $self->{'collection'}{$name} }, $value;
            }
        }
        $self->{'collection'}{$name} =
            [ sort { $a <=> $b } @{ $self->{'collection'}{$name} } ];
    }
    return;
}

sub remove
{
    my ( $self, $bond_parameters ) = @_;
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
