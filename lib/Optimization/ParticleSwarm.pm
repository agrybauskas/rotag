package ParticleSwarm;

use strict;
use warnings;

use Optimization::Particle;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $parameters, $particle_num, $options ) = @_;
    my ( $seed ) = $options->{'seed'};

    $seed //= 23;

    srand( $seed );

    my $self = { 'particles' => undef,
                 'cost_function' => undef,
                 'optimal_value' => undef,
                 'optimal_parameters' => undef };
    for my $i ( 0..$particle_num-1 ) {
        my $id = $i + 1;
        my $particle = Particle->new( $parameters );
        if( $particle->no_values ) {
            for my $name ( keys %{ $particle->{'parameters'} } ) {
                my $parameter = $particle->{'parameters'}{$name};
                my $min = $parameter->min;
                my $max = $parameter->max;
                $parameter->value( $min + rand( $max - $min ) );
            }
        }
        $self->{'particles'}{$id} = $particle;
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_cost_function
{
    my ( $self, $cost_function ) = @_;
    $self->{'cost_function'} = $cost_function;
}

# --------------------------------- Methods ----------------------------------- #

sub optimize
{
    my ( $self, $iterations ) = @_;

    my $cost_function = $self->{'cost_function'};
    if( ! defined $cost_function ) {
        die "'cost_function' value for Optimization::ParticleSwarm is " .
            "mandatory.\n";
    }

    my $particles = $self->{'particles'};
    for my $i ( 0..$iterations-1 ) {
        for my $id ( keys %{ $particles } ) {
            my $particle = $particles->{$id};
            my $parameters = $particle->{'parameters'};
            my $speed = $particle->speed;

            if( defined $speed ) {
            }

            $particle->value( $cost_function->( $parameters ) );

            if( ! defined $self->{'optimal_value'} ||
                $particle->value < $self->{'optimal_value'} ) {
                $self->{'optimal_value'} = $particle->value;
                $self->{'optimal_parameters'} = $parameters;
            }
        }

        # Sets speed for next iteration.
        for my $id ( keys %{ $particles } ) {
            my $particle = $particles->{$id};
            my $parameters = $particle->{'parameters'};

            my $weight = $particle->value / $self->{'optimal_value'};

            for my $name ( keys %{ $parameters } ) {
                my $parameter = $parameters->{$name};
                my $optimal_parameter = $self->{'optimal_parameters'}{$name};
                # $parameter->speed($optimal_parameter->value-$parameter->value);
            }
        }
    }
}

1;
