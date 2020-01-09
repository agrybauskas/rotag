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
                 'global_optimal_value' => undef,
                 'global_optimal_param' => undef,
                 'local_optimal_value' => undef,
                 'local_optimal_param' => undef };
    for my $i ( 0..$particle_num-1 ) {
        my $id = $i + 1;
        my $particle = Particle->new( $parameters );
        if( $particle->no_values ) {
            for my $name ( sort keys %{ $particle->{'parameters'} } ) {
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

sub global_optimal_value
{
    my ( $self, $global_optimal_value ) = @_;
    if( scalar @_ == 2 ) {
        $self->{'global_optimal_value'} = $global_optimal_value;
    }
    return $self->{'global_optimal_value'};
}

sub global_optimal_param
{
    my ( $self, $global_optimal_param ) = @_;
    if( scalar @_ == 2 ) {
        $self->{'global_optimal_param'} = $global_optimal_param;
    }
    return $self->{'global_optimal_param'};
}

sub local_optimal_value
{
    my ( $self, $local_optimal_value ) = @_;
    if( scalar @_ == 2 ) {
        $self->{'local_optimal_value'} = $local_optimal_value;
    }
    return $self->{'local_optimal_value'};
}

sub local_optimal_param
{
    my ( $self, $local_optimal_param ) = @_;
    if( scalar @_ == 2 ) {
        $self->{'local_optimal_param'} = $local_optimal_param;
    }
    return $self->{'local_optimal_param'};
}

# --------------------------------- Methods ----------------------------------- #

sub optimize
{
    my ( $self, $iterations, $options ) = @_;

    my $cost_function = $self->{'cost_function'};
    if( ! defined $cost_function ) {
        die "'cost_function' value for Optimization::ParticleSwarm is " .
            "mandatory.\n";
    }

    my $particles = $self->{'particles'};
    for my $i ( 0..$iterations-1 ) {
        for my $id ( sort keys %{ $particles } ) {
            my $particle = $particles->{$id};
            my $parameters = $particle->{'parameters'};

            if( ! defined $particle->speed ) {
                my %speed = ();
                for my $key ( sort keys %{ $parameters } ) {
                    $speed{$key} =
                        rand( 1 ) *
                        ( $parameters->{$key}->uniform() -
                          $parameters->{$key}->value );
                }
                $particle->speed( \%speed );
            }

            for my $key ( sort keys %{ $parameters } ) {
                my $parameter_value = $parameters->{$key}->value;
                $parameters->{$key}->value(
                    $parameter_value + $particle->speed->{$key}
                );
            }

            $particle->value( $cost_function->( $parameters ) );

            if( ! defined $self->{'global_optimal_value'} ||
                $particle->value <= $self->{'global_optimal_value'} ) {
                $self->{'global_optimal_value'} = $particle->value;
                $self->{'global_optimal_param'} = $parameters;
            }
            if( ! exists $self->{'local_optimal_value'}{$id} ||
                $particle->value <= $self->{'local_optimal_value'}{$id} ) {
                $self->{'local_optimal_value'}{$id} = $particle->value;
                $self->{'local_optimal_param'}{$id} = $parameters;
            }
        }

        # Sets speed for next iteration.
        for my $id ( sort keys %{ $particles } ) {
            my $particle = $particles->{$id};
            my $parameters = $particle->{'parameters'};
            my %updated_speed = ();
            for my $key ( sort keys %{ $parameters } ) {
                my $parameter = $parameters->{$key};
                my $global_optimal_parameter =
                    $self->{'global_optimal_param'}{$key};
                my $local_optimal_parameter =
                    $self->{'local_optimal_param'}{$id}{$key};
                $updated_speed{$key} =
                    rand( 1 ) *
                    ( ( $global_optimal_parameter->value - $parameter->value ) -
                      ( $local_optimal_parameter->value - $parameter->value ) );
            }
            $particle->speed( \%updated_speed );
        }
    }
}

1;
