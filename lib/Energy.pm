package Energy;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my $self = { 'energy_type' => undef,
                 'atoms' => undef,
                 'value' => undef };
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_energy
{
    my ( $self, $energy_type, $atoms, $value ) = @_;
    $self->{'energy_type'} = $energy_type;
    $self->{'atoms'} = $atoms;
    $self->{'value'} = $value;
}

sub energy_type
{
    my ( $self ) = @_;
    return $self->{'energy_type'};
}

sub atoms
{
    my ( $self ) = @_;
    return $self->{'atoms'};
}

sub value
{
    my ( $self ) = @_;
    return $self->{'value'};
}

# --------------------------------- Methods ----------------------------------- #

1;
