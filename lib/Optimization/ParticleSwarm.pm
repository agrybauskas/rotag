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

    my $self = { 'particles' => undef };
    for my $i ( 0..$particle_num-1 ) {
        my $id = $i + 1;
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

1;
