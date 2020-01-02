package Optimization;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( particle_swarm );

use Optimization::Parameters;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $particles ) = @_;
    my $self = {
        'particles' => undef,
        'cost_function' => undef,
    };

    for my $name ( keys %{ $particles } ) {
        my $particle = Parameters->new( {
            'key' => $name,
            'min_range' => $particles->{$name}{'min_range'},
            'max_range' => $particles->{$name}{'max_range'}
        } );
        $self->{'particles'}{$name} = $particle;
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_cost_function
{
    my ( $self, $cost_function ) = @_;
    $self->{'cost_function'} = $cost_function;
}

# ------------------------- Optimization algorithms --------------------------- #

sub particle_swarm
{
    my ( $particles, $options ) = @_;
    my $cost_function = $particles->{'cost_function'};
    my ( $seed ) = ( $options->{'seed'} );
    $seed //= 23;

    srand( $seed );

    if( ! defined $cost_function ) {
        die "Cost function is missing. It has to be set.\n";
    }

    for my $key ( keys %{ $particles->{'particles'} } ) {
        my ( $particle ) = $particles->{'particles'}{$key};

        if( ! defined $particle->value ) {
            $particle->value(
                $particle->min_range +
                rand( $particle->max_range - $particle->min_range )
            );
        }
    }
}

1;
