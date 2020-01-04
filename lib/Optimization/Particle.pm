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

1;
