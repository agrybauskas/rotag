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
    for my $bond_parameter ( @{ $bond_parameters } ) {
        my ( $name ) = keys %{ $bond_parameter };
        my @values = @{ $bond_parameter->{$name} };
        while( @values ) {
            my ( $value ) = shift @values;
            my @del_indexes =
                grep { $self->{'collection'}{$name}[$_] eq $value }
                     ( 0..$#{ $self->{'collection'}{$name} } );
            for my $del_index ( @del_indexes ) {
                splice @{ $self->{'collection'}{$name} }, $del_index, 1;
            }
        }
    }
    return;
}

# --------------------------------- Methods ----------------------------------- #

sub get_names
{
    my ( $self ) = @_;
    return $self->{'order'};
}

sub get_values
{
    my ( $self, $names, $options ) = @_;
    my $cache = ( $options->{'cache'} );
    $names //= $self->{'order'};
    $cache //= 0;
    my $permuted_values = permutation(
        scalar( @{ $names } ),
        [],
        [ map { $self->{'collection'}{$_} } @{ $names } ],
        []
        );
    if( $cache ) {
        my $names_key = join ',', @{ $names };
        $self->{'cache'}{$names_key} = $permuted_values;
    }
    return $permuted_values;
}

sub get_cache
{
    my ( $self, $names ) = @_;
    $names //= $self->{'order'};
    my $names_key = join ',', @{ $names };
    return $self->{'cache'}{$names_key} if exists $self->{'cache'}{$names_key};
    return [];
}

sub exists
{
    my ( $self, $name ) = @_;
    return 1 if defined $name && exists $self->{'collection'}{$name};
    return 0;
}

1;
