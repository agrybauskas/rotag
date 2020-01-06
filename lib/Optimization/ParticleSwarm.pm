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
                 'optimal_position' => undef,
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
        # TODO: should be added more than min function when searching for
        # optimal solution.
        my $best_position;
        my $best_value;
        for my $id ( keys %{ $particles } ) {
            my $particle = $particles->{$id};
            my $parameters = $particle->{'parameters'};
            my @parameter_values =
                map { $parameters->{$_}->value } sort keys %{ $parameters };
            my $best_particle_value = $cost_function->( @parameter_values );

            $particle->position( $best_particle_value );

            # Check for best overal value in all particles.
            if( ( ! defined $best_position || ! defined $best_value ) ||
                ( $best_particle_value <= $best_value  ) ) {
                $best_position = \@parameter_values;
                $best_value = $best_particle_value;
            }
        }
    }
}

1;
