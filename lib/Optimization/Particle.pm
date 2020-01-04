package Particle;

use strict;
use warnings;

use Optimization::Parameter;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $parameters, $options ) = @_;

    my $self = { 'parameters' => undef,
                 'position' => undef,
                 'speed' => undef };

    for my $name ( keys %{ $parameters } ) {
        my $parameter = Parameter->new( {
            'key' => $name,
            'min' => $parameters->{$name}{'min'},
            'max' => $parameters->{$name}{'max'},
            'value' => undef,
        } );
        $self->{'parameters'}{$name} = $parameter;
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub get_values
{
    my ( $self ) = @_;
    my %values = ();
    for my $name ( keys %{ $self->{'parameters'} } ) {
        $values{$name} = $self->{'parameters'}{$name}->value;
    }
    return \%values;
}

sub position
{
    my ( $self, $position ) = @_;
    if( scalar @_ == 2  ) {
        $self->{'position'} = $position;
    }
    return $self->{'position'};
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

sub no_values
{
    my ( $self ) = @_;
    my $parameter_values = $self->get_values;
    if( map { $parameter_values->{$_} } keys %{ $parameter_values } ) {
        return 1;
    }
    return 0;
}

1;
