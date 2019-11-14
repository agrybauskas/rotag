package PredictSidechains;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sidechain_positions );

use Graph;
use GraphViz;

use PDBxParser qw( filter_new
                   unique_residue_key );
use Grid qw( grid_box
             identify_neighbour_cells );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub sidechain_positions
{
    my ( $args ) = @_;
    my ( $parameters, $atom_site ) =
        ( $args->{'parameters'}, $args->{'atom_site'} );

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
    # my $graph_viz = GraphViz->new(); # NOTE: only for development purposes.

    for my $cell ( keys %{ $grid_box_cas } ) {
        my $neighbour_cell_atom_ids = $neighbouring_cells->{$cell};
        for my $atom_id ( @{ $grid_box_cas->{$cell} } ) {
            my $unique_residue_key =
                unique_residue_key( $atom_site_cas->{$atom_id} );

            $interaction_graph->add_vertex( $unique_residue_key );
            # $graph_viz->add_node( $unique_residue_key );

            my $neighbour_atom_ids =
                [ grep { $atom_id ne $_ } @{ $grid_box_cas->{$cell} }  ];
            push @{ $neighbour_atom_ids }, @{ $neighbour_cell_atom_ids };

            for my $neighbour_atom_id ( @{ $neighbour_atom_ids } ) {
                my $neighbour_residue_key =
                    unique_residue_key( $atom_site_cas->{$neighbour_atom_id} );

                $interaction_graph->add_edge( $unique_residue_key,
                                              $neighbour_residue_key );
                # $graph_viz->add_node( $neighbour_residue_key );
                # $graph_viz->add_edge( $unique_residue_key =>
                #                           $neighbour_residue_key );
            }
        }
    }

    use Data::Dumper;
    print Dumper $interaction_graph;
    # print $graph_viz->as_png;
}

# ----------------------------------------------------------------------------- #

1;
