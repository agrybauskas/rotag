package Parameter;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my $self = {};

    # Mandatory arguments.
    $self->{'key'} = $args->{'key'};
    $self->{'min_range'} = $args->{'min_range'};
    $self->{'max_range'} = $args->{'max_range'};

    if( ! defined $self->{'key'} ) {
        die "'key' value for Optimization::Parameter is mandatory.\n";
    }
    if( ! defined $self->{'min_range'} ) {
        die "'min_range' value for Optimization::Parameter is mandatory.\n";
    }
    if( ! defined $self->{'max_range'} ) {
        die "'max_range' value for Optimization::Parameter is mandatory.\n";
    }

    # Optional arguments.
    $self->{'value'} = $args->{'value'};
    $self->{'speed'} = $args->{'speed'};

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub key
{
    my ( $self ) = @_;
    return $self->{'key'};
}

sub min_range
{
    my ( $self ) = @_;
    return $self->{'min_range'};
}

sub max_range
{
    my ( $self ) = @_;
    return $self->{'max_range'};
}

sub value
{
    my ( $self, $value ) = @_;
    if( scalar @_ == 2  ) {
        $self->{'value'} = $value;
    }
    return $self->{'value'};
}

sub speed
{
    my ( $self, $speed ) = @_;
    if( scalar @_ == 2  ) {
        $self->{'speed'} = $speed;
    }
    return $self->{'speed'};
}

# --------------------------------- Methods ----------------------------------- #

1;
