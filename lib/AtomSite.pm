package AtomSite;

use strict;
use warnings;

use Grid qw( grid_box );

# ----------------------- Constructors and Destructors ------------------------ #

sub new
{
    my ( $class, $args )  = @_;

    my %self = (
        # { <id> => { 'id' => <id> , 'Cartn_x' => <Cartn_x>, ... }, ... }
        'atoms'          => $args->{'atoms'},
        # { <atoms_id> => [ <atoms_id>, ... ], ... }
        'connections'    => $args->{'connections'},
        # { <atom_id> => <hybridization>, ... }
        'hybridizations' => $args->{'interactions'},
        # { <atom_id> => { <atom_id> => { <type> => <value> }, ... }, ... }
        'interactions'   => $args->{'interactions'},
        # { "<cell_x_id>,<cell_y_id>,<cell_z_ud>" => [ <atom_id>, ... ], ... }
        'grid_box'       => $args->{'grid_box'}
    );

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
