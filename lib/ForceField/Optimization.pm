package Optimization;

use strict;
use warnings;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;

    my $self = { 'parameters'  => undef,
                 'calculation' => undef,
                 'estimation'  => undef,
                 'constraints' => undef };

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub parameters
{

}

# --------------------------------- Methods ----------------------------------- #

sub calculation
{

}

sub estimation
{

}

sub contraints
{

}

1;
