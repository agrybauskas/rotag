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
                 'cost_function' => undef };
    for my $i ( 0..$particle_num-1 ) {
        my $id = $i + 1;
        my $particle = Particle->new( $parameters );
        if( $particle->no_values ) {
            for my $name ( keys %{ $particle->{'parameters'} } ) {
                my $parameter = $particle->{'parameters'}{$name};
                my $min = $parameter->min;
                my $max = $parameter->max;
                $parameter->value( $min + rand( $min - $max ) );
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

    if( ! defined $self->{'cost_function'} ) {
        die "'cost_function' value for Optimization::ParticleSwarm is " .
            "mandatory.\n";
    }

    for my $i ( 0..$iterations-1 ) {
        print $i, "\n";
    }
}

1;
