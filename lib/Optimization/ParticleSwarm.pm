package ParticleSwarm;

use strict;
use warnings;

use Optimization::Particle;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $parameters, $particle_num, $options ) = @_;

    my $self = {};

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

1;
