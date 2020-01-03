package Optimization;

use strict;
use warnings;

use Optimization::Parameter;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $parameters, $particle_num, $options ) = @_;
    my ( $seed ) = ( $options->{'seed'} );
    $seed //= 23;

    my $self = {
        'particles' => undef,
        'cost_function' => undef,
    };

    srand( $seed );

    for my $i ( 0..$particle_num-1 ) {
        my $id = $i+1;
        for my $name ( keys %{ $parameters } ) {
            my $parameter = Parameters->new( {
                'key' => $name,
                'min_range' => $parameters->{$name}{'min_range'},
                'max_range' => $parameters->{$name}{'max_range'},
                'value' =>
                    $parameters->{$name}{'min_range'} +
                    rand( $parameters->{$name}{'max_range'} -
                          $parameters->{$name}{'min_range'} )
            } );
            $self->{'particles'}{$id}{'parameters'}{$name} = $parameter;
        }
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_cost_function
{
    my ( $self, $cost_function ) = @_;
    $self->{'cost_function'} = $cost_function;
}

1;
