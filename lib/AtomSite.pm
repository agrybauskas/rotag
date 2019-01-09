package AtomSite;

use strict;
use warnings;

use Grid qw( grid_box );

# ----------------------- Constructors and Destructors ------------------------ #

sub new
{
    my ( $class, $args )  = @_;

    my %self = ( 'atoms'        => $args->{'atoms'},
                 'connections'  => $args->{'connections'},
                 'interactions' => $args->{'interactions'},
                 'grid_box'     => $args->{'grid_box'} );

    return bless \%self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_grid_box
{
    my ( $self, $edge_length ) = @_;

    ( $self->{'grid_box'} ) = grid_box( $self->{'atoms'}, $edge_length );

    return;
}

sub get_grid_box
{
    my ( $self ) = @_;

    return $self->{'grid_box'};
}

# --------------------------------- Methods ----------------------------------- #

# ----------------------------- Static functions ------------------------------ #

1;
