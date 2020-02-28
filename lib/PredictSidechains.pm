package PredictSidechains;

use strict;
use warnings;

use Graph;

use PDBxParser qw( filter_new
                   unique_residue_key );
use Grid qw( grid_box
             identify_neighbour_cells );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my $self = {
        'atom_site' => $args->{'atom_site'},
        'rotamer_angles' => $args->{'rotamer_angles'},
        'rotamer_energies' => $args->{'rotamer_energies'},
        'interaction_graph' => undef,
        'parameters' => $args->{'parameters'}
    };
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

# --------------------------------- Methods ----------------------------------- #

sub interaction_graph
{
    my ( $self ) = @_;

    my ( $parameters, $atom_site ) =
        ( $self->{'parameters'}, $self->{'atom_site'} );

    if( ! defined $parameters ) {
        die "parameters are not set.\n";
    }
    if( ! defined $atom_site ) {
        die "atom site is not supplied.\n";
    }

    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};

    # Chooses only CA atoms, because from them, boundary interaction conditions
    # are measured.
    my $atom_site_cas = filter_new( $atom_site,
                                    { 'include' =>
                                          { 'label_atom_id' => [ 'CA' ] } } );

    my ( $grid_box_cas ) = grid_box( $parameters, $atom_site_cas,
                                     $edge_length_interaction );
    my $neighbouring_cells = identify_neighbour_cells( $grid_box_cas );

    my $interaction_graph = Graph->new();
    for my $cell ( keys %{ $grid_box_cas } ) {
        my $neighbour_cell_atom_ids = $neighbouring_cells->{$cell};
        for my $atom_id ( @{ $grid_box_cas->{$cell} } ) {
            my $unique_residue_key =
                unique_residue_key( $atom_site_cas->{$atom_id} );

            $interaction_graph->add_vertex( $unique_residue_key );

            my $neighbour_atom_ids =
                [ grep { $atom_id ne $_ } @{ $neighbour_cell_atom_ids }  ];

            for my $neighbour_atom_id ( @{ $neighbour_atom_ids } ) {
                my $neighbour_residue_key =
                    unique_residue_key( $atom_site_cas->{$neighbour_atom_id} );

                $interaction_graph->add_vertex( $unique_residue_key );
                $interaction_graph->add_edge( $unique_residue_key,
                                              $neighbour_residue_key );
            }
        }
    }

    $self->{'interaction_graph'} = $interaction_graph;

    return;
}

sub choose
{
    my ( $self ) = @_;

    my ( $rotamer_angles, $rotamer_energies, $interaction_graph ) = (
        $self->{'rotamer_angles'},
        $self->{'rotamer_energies'},
        $self->{'interaction_graph'}
    );

    if( ! defined $rotamer_angles ) {
        die "rotamer angles are not supplied by '_[local]_rotamer_angle' tag.\n";
    }
    if( ! defined $rotamer_energies ) {
        die "rotamer energies are not supplied by '_[local]_rotamer_energy'" .
            " tag.\n";
    }
    if( ! defined $interaction_graph ) {
        $self->interaction_graph();
    }

    # Generates related data structure that joins rotamer angle with energy data.
    # TODO: should use similar function to related_category_data().
    for my $rotamer_angle_id ( keys %{ $rotamer_angles } ) {
        my $rotamer_angle = $rotamer_angles->{$rotamer_angle_id};
        my $rotamer_angle_type = $rotamer_angle->{'type'};
        my $rotamer_id = $rotamer_angle->{'rotamer_id'};
        $rotamer_energies->{$rotamer_id}{'angles'}{$rotamer_angle_type} =
            $rotamer_angle;
        # TODO: removing frequencies, because they can be changed and, for now,
        # they are ignored.
        delete $rotamer_energies->{$rotamer_id}{'angles'}{$rotamer_angle_type}
                                                         {'frequency'};
    }
}

1;
