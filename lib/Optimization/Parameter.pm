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
    $self->{'min'} = $args->{'min'};
    $self->{'max'} = $args->{'max'};

    if( ! defined $self->{'key'} ) {
        die "'key' value for Optimization::Parameter is mandatory.\n";
    }
    if( ! defined $self->{'min'} ) {
        die "'min' value for Optimization::Parameter is mandatory.\n";
    }
    if( ! defined $self->{'max'} ) {
        die "'max' value for Optimization::Parameter is mandatory.\n";
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub key
{
    my ( $self ) = @_;
    return $self->{'key'};
}

sub value
{
    my ( $self, $value ) = @_;
    if( scalar @_ == 2 ) {
        $self->{'value'} = $value;
    }
    return $self->{'value'};
}

sub min
{
    my ( $self ) = @_;
    return $self->{'min'};
}

sub max
{
    my ( $self ) = @_;
    return $self->{'max'};
}

sub uniform
{
    my ( $parameter ) = @_;
    my $min = $parameter->min;
    my $max = $parameter->max;
    return $min + rand( $max - $min );
}

1;
