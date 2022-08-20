package PredictSidechains;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( predict_sidechains );

use Clone qw( clone );
use Graph::Undirected;

use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use ForceField::Parameters;
use ForceField::Bonded qw( general );
use ForceField::NonBonded qw( general );
use Measure qw( bond_length
                energy );
use PDBxParser qw( extract
                   filter_new
                   filter_by_unique_residue_key
                   unique_residue_key
                   to_pdbx );
use PseudoAtoms qw( pairwise_rotamer_energy
                    replace_with_rotamer );
use Grid qw( grid_box
             identify_neighbour_cells );
use SidechainModels qw( rotation_only );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------- Simple functions ------------------------------ #

sub predict_sidechains
{
    my ( $args ) = @_;

    my ( $parameters, $atom_site, $rotamer_energies, $rotamer_angles,
         $non_bonded_potential, $bonded_potential, $options ) = (
        $args->{'parameters'}, $args->{'atom_site'}, $args->{'rotamer_energies'},
        $args->{'rotamer_angles'}, $args->{'non_bonded_potential'},
        $args->{'bonded_potential'}, $args->{'options'}
    );

    # Error messages for missing arguments.
    if( ! defined $atom_site ) {
        die "atom site is not supplied.\n";
    }
    if( ! defined $rotamer_angles ) {
        die "rotamer angles are not supplied by '_[local]_rotamer_angle' tag.\n";
    }
    if( ! defined $rotamer_energies ) {
        die "rotamer energies are not supplied by '_[local]_rotamer_energy'" .
            " tag.\n";
    }
    if( ! defined $parameters ) {
        die "parameters are not set.\n";
    }

    # Generates look up table for easier data accesibility.
    my %rotamer_to_residue = ();
    my %residue_to_rotamer = ();
    my %rotamer_to_angles = ();
    for my $rotamer_angle_id ( keys %{ $rotamer_angles } ) {
        my $rotamer_angle = $rotamer_angles->{$rotamer_angle_id};
        my $rotamer_id = $rotamer_angle->{'rotamer_id'};
        my $unique_residue_key =
            sprintf '%s,%s,%s,%s',
            $rotamer_angle->{'label_seq_id'},
            $rotamer_angle->{'label_asym_id'},
            $rotamer_angle->{'pdbx_PDB_model_num'},
            $rotamer_angle->{'label_alt_id'};

        $rotamer_to_residue{$rotamer_id}{$unique_residue_key} = 1;
        $residue_to_rotamer{$unique_residue_key}{$rotamer_id} = 1;
        $rotamer_to_angles{$rotamer_id}{$rotamer_angle_id} = 1;
    }

    # Determining interaction grid.
    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};
    my $cutoff_atom =
        $parameters->{'_[local]_force_field'}{'cutoff_atom'};

    # Chooses only CA atoms, because from them, boundary interaction conditions
    # are measured.
    my $atom_site_cas = filter_new( $atom_site,
                                    { 'include' =>
                                          { 'label_atom_id' => [ 'CA' ] } } );
    my ( $grid_box_cas, $grid_ca_atom_pos ) =
        grid_box( $parameters, $atom_site_cas, $edge_length_interaction,
                  extract( $atom_site_cas,
                           { 'data' => [ 'id' ], 'is_list' => 1 } ) );
    my $neighbouring_cells_cas =
        identify_neighbour_cells( $grid_box_cas, $grid_ca_atom_pos );

    my %residue_pairs = ();
    my %residue_to_grid = ();
    for my $grid_id ( keys %{ $grid_ca_atom_pos } ) {
        for my $atom_id ( @{ $grid_ca_atom_pos->{$grid_id} } ) {
            for my $neighbour_atom_id ( @{ $neighbouring_cells_cas->{$grid_id} } ) {
                next if $atom_id eq $neighbour_atom_id ||
                    ( defined $residue_pairs{$atom_id} &&
                      defined $residue_pairs{$atom_id}
                                            {$neighbour_atom_id} &&
                      $residue_pairs{$atom_id}
                                    {$neighbour_atom_id} );

                next if bond_length(
                    [ [ map { $atom_site->{$atom_id}{$_} }
                          ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ],
                      [ map { $atom_site->{$neighbour_atom_id}{$_} }
                          ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ] ] ) >=
                    $edge_length_interaction;

                $residue_pairs{$atom_id}
                              {$neighbour_atom_id} = 1;
            }
            $residue_to_grid{$atom_id} = $grid_id;
        }
    }

    # Generates graph for protein.
    # TODO: make sure that interactions are in both ways.
    my $interaction_graph = new Graph::Undirected;
    for my $ca_atom_id ( keys %residue_pairs ) {
        if( ! $interaction_graph->has_vertex( $ca_atom_id ) ) {
            $interaction_graph->add_vertex( $ca_atom_id );
            $interaction_graph->set_vertex_attribute( $ca_atom_id, 'visited', 0);
        }

        for my $neighbour_ca_atom_id (
            keys %{ $residue_pairs{$ca_atom_id} } ) {
            if( ! $interaction_graph->has_vertex( $neighbour_ca_atom_id ) ) {
                $interaction_graph->add_vertex( $neighbour_ca_atom_id );
                $interaction_graph->set_vertex_attribute( $neighbour_ca_atom_id,
                                                          'visited', 0);
            }

            if( ! $interaction_graph->has_edge( $ca_atom_id,
                                                $neighbour_ca_atom_id ) ) {
                $interaction_graph->add_edge( $ca_atom_id,
                                              $neighbour_ca_atom_id );
                $interaction_graph->set_edge_attribute( $ca_atom_id,
                                                        $neighbour_ca_atom_id,
                                                        'visited', 0);
            }
        }
    }

    # Starts at random vertex and travels bread-first.
    # TODO: there will definately be problems with unconnected graphs. Have to
    # enter the second graph and etc.
    my $start_node = $interaction_graph->random_vertex;
    my @next_nodes = $interaction_graph->neighbours( $start_node );

    # Keeps structure data if it was already calculated.
    my %rotamer_to_atom_site = ();
    my %rotamer_interaction_counter = (); # Will be used to remove rotamers with
                                          # no interacting candidates.
    my %predicted_rotamer_pairs = ();

    while( @next_nodes ) {
        for my $node ( @next_nodes ) {
            my $unique_residue_key = unique_residue_key($atom_site_cas->{$node});
            my @rotamer_ids =
                keys %{ $residue_to_rotamer{$unique_residue_key} };
            my $rotamer_site =
                filter_by_unique_residue_key( $atom_site,
                                              $unique_residue_key,
                                              1 );

            connect_atoms( $parameters, $rotamer_site );
            hybridization( $parameters, $rotamer_site );
            rotation_only( $parameters, $rotamer_site );

            for my $rotamer_id ( @rotamer_ids ) {
                my @angle_ids = keys %{ $rotamer_to_angles{$rotamer_id} };
                my %angles =
                    map { $rotamer_angles->{$_}{'type'} =>
                          $rotamer_angles->{$_}{'value'} }
                        @angle_ids;

                if( ! defined $rotamer_interaction_counter{$rotamer_id} ) {
                    $rotamer_interaction_counter{$rotamer_id} = 0;
                }

                my %rotamer_site;
                if( ! defined $rotamer_to_atom_site{$rotamer_id} ) {
                    %rotamer_site = %{ clone( $rotamer_site ) };
                    replace_with_rotamer( $parameters, \%rotamer_site,
                                          $unique_residue_key, \%angles );
                    $rotamer_to_atom_site{$rotamer_id} = { %rotamer_site };
                } else {
                    %rotamer_site = %{ $rotamer_to_atom_site{$rotamer_id} };
                }

                for my $neighbour_node ( keys %{ $residue_pairs{$node} } ) {
                    my $neighbour_unique_residue_key =
                        unique_residue_key( $atom_site_cas->{$neighbour_node} );
                    my @neighbour_rotamer_ids =
                        keys %{$residue_to_rotamer{$neighbour_unique_residue_key}};
                    my $neighbour_rotamer_site = filter_by_unique_residue_key(
                        $atom_site,
                        $neighbour_unique_residue_key,
                        1
                    );

                    # TODO: move code so, it would be calculated once.
                    connect_atoms( $parameters, $neighbour_rotamer_site );
                    hybridization( $parameters, $neighbour_rotamer_site );
                    rotation_only( $parameters, $neighbour_rotamer_site );

                    for my $neighbour_rotamer_id ( @neighbour_rotamer_ids ) {
                        my @neighbour_angle_ids =
                            keys %{ $rotamer_to_angles{$neighbour_rotamer_id} };
                        my %neighbour_angles =
                            map { $rotamer_angles->{$_}{'type'} =>
                                  $rotamer_angles->{$_}{'value'} }
                                @neighbour_angle_ids;

                        my %neighbour_rotamer_site;
                        if(!defined $rotamer_to_atom_site{$neighbour_rotamer_id}){
                            %neighbour_rotamer_site =
                                %{ clone( $neighbour_rotamer_site ) };
                            replace_with_rotamer( $parameters,
                                                  \%neighbour_rotamer_site,
                                                  $neighbour_unique_residue_key,
                                                  \%angles );
                            $rotamer_to_atom_site{$neighbour_rotamer_id} =
                                { %neighbour_rotamer_site };
                        } else {
                            %neighbour_rotamer_site =
                                %{ $rotamer_to_atom_site{$neighbour_rotamer_id}};
                        }

                        # Calculate pairwise energy.
                        my $pairwise_energy_sum = pairwise_rotamer_energy(
                            $parameters,
                            \%rotamer_site,
                            \%neighbour_rotamer_site,
                            \&ForceField::Bonded::general,
                            \&ForceField::NonBonded::general,
                        );

                        # Does not reach cut off limit.
                        if( $pairwise_energy_sum <= $cutoff_atom  ) {
                            # Stores energy value.
                            $predicted_rotamer_pairs{$rotamer_id}
                                                    {$neighbour_rotamer_id} =
                                $pairwise_energy_sum;
                            $predicted_rotamer_pairs{$neighbour_rotamer_id}
                                                    {$rotamer_id} =
                                $pairwise_energy_sum;

                            # Increases the counter in order to track removable
                            # rotamers.
                            $rotamer_interaction_counter{$rotamer_id} += 1;
                            $rotamer_interaction_counter{$neighbour_rotamer_id}
                                                                            += 1;
                        }
                    }
                }
            }
        }

        @next_nodes = (); # Temporare reset.
    }

    return \%predicted_rotamer_pairs, $interaction_graph;
}

1;
