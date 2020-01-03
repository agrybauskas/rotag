package Optimization;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( particle_swarm );

use Optimization::Parameter;

# ------------------------- Optimization algorithms --------------------------- #

sub particle_swarm
{
    my ( $particles, $options ) = @_;
    my $cost_function = $particles->{'cost_function'};

    if( ! defined $cost_function ) {
        die "Cost function is missing. It has to be set.\n";
    }

    for my $id ( keys %{ $particles->{'particles'} } ) {
        my $particle = $particles->{'particles'}{$id};
    }
}

1;
