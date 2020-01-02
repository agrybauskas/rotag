package Optimization;

use strict;
use warnings;

use Optimization::Parameters;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $particles ) = @_;
    my $self = {
        'particles' => undef,
        'min_value' => undef
    };

    for my $name ( keys %{ $particles } ) {
        my $particle = Parameters->new( {
            'key' => $name,
            'min_range' => $particles->{$name}{'min_range'},
            'max_range' => $particles->{$name}{'max_range'}
        } );
        $self->{'particles'} = $particle;
    }

    return bless $self, $class;
}

# ------------------------- Optimization algorithms --------------------------- #

sub particle_swarm
{
    my ( $particles ) = @_;
}


1;
