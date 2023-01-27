package BondParameter;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class,  ) = @_;

    my $self = {};

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub get_parameter_atom_ids
{
    my ( $self, $atom_id, $parameter_name ) = @_;
    return $self->{'id'}{$atom_id}{$parameter_name}{'atom_ids'};
}

sub get_ordered_parameter_names
{
    my ( $self, $atom_id ) = @_;
    return [
        sort { $self->{'id'}{$atom_id}{$a}{'order'} <=>
               $self->{'id'}{$atom_id}{$b}{'order'} }
        keys %{ $self->{'id'}{$atom_id} }
    ];
}

sub get_parameters
{
    my ( $self, $residue_unique_key ) = @_;
    return {} if ! exists $self->{'residue_unique_key'}{$residue_unique_key};
    return $self->{'residue_unique_key'}{$residue_unique_key};
}

sub get_parameter
{
    my ( $self, $residue_unique_key, $parameter_name ) = @_;
    return $self->{'residue_unique_key'}{$residue_unique_key}{$parameter_name};
}

sub get_residue_unique_keys
{
    my ( $self ) = @_;
    return [] if ! exists $self->{'residue_unique_key'};
    return [ keys %{ $self->{'residue_unique_key'} } ];
}

1;
