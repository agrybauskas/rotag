package AtomSite;

use strict;
use warnings;

use ConnectAtoms qw( connect_atoms_new );
use Grid qw( grid_box
             identify_neighbour_cells_new );

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
        'hybridizations' => $args->{'hybridizations'},
        # { <atom_id> => { <atom_id> => { <type> => <value> }, ... }, ... }
        'interactions'   => $args->{'interactions'},
        # { "<cell_x_id>,<cell_y_id>,<cell_z_id>" =>
        #       { 'atom_ids' => [ <atom_id>, ... ],
        #         'neighbours' => [ "<cell_x_id>,<cell_y_id>,<cell_z_id>", ... ] } }
        'grid_box'       => $args->{'grid_box'},
    );

    return bless \%self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub set_grid_box
{
    my ( $self, $options ) = @_;
    my $edge_length = $options->{'edge_length'};

    my ( $grid_box ) = grid_box( $self->{'atoms'}, $edge_length );
    my $neighbouring_cells = identify_neighbour_cells_new( $grid_box );

    for my $cell ( keys %{ $grid_box } ) {
        $self->{'grid_box'}{$cell}{'atom_ids'} = $grid_box->{$cell};
        $self->{'grid_box'}{$cell}{'neighbours'} = $neighbouring_cells->{$cell};
    }

    return;
}

sub get_grid_box
{
    my ( $self ) = @_;

    return $self->{'grid_box'};
}

sub set_connections
{
    my ( $self, $options ) = @_;
    my $edge_length = $options->{'edge_length'};

    if( ! defined $self->{'grid_box'} ) {
        $self->set_grid_box( { 'edge_length' => $edge_length } );
    }

    $self->{'connections'} =
        connect_atoms_new( $self->{'atoms'}, $self->{'grid_box'} );

    return;
}

sub get_connections
{
    my ( $self ) = @_;

    return $self->{'connections'};
}

# --------------------------------- Methods ----------------------------------- #

# ----------------------------- Static functions ------------------------------ #

1;
