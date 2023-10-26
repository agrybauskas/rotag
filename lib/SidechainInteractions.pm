package SidechainInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( predict_sidechains );

use Clone qw( clone );

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

sub new
{
    my ( $class, $parameters, $atom_site, $rotamer_angles,
         $rotamer_energies ) = @_;

    # if( ! defined $atom_site ) {
    #     die "atom site is not supplied.\n";
    # }
    # if( ! defined $rotamer_angles ) {
    #     die "rotamer angles are not supplied by '_[local]_rotamer_angle' tag.\n";
    # }
    # if( ! defined $rotamer_energies ) {
    #     die "rotamer energies are not supplied by '_[local]_rotamer_energy'" .
    #         " tag.\n";
    # }

    my $self = {};

    # my %rotamer_to_angles = ();
    # for my $rotamer_angle_id ( keys %{ $rotamer_angles } ) {
    #     my $rotamer_angle = $rotamer_angles->{$rotamer_angle_id};
    #     my $rotamer_id = $rotamer_angle->{'rotamer_id'};
    #     my $unique_residue_key =
    #         sprintf '%s,%s,%s,%s',
    #         $rotamer_angle->{'label_seq_id'},
    #         $rotamer_angle->{'label_asym_id'},
    #         $rotamer_angle->{'pdbx_PDB_model_num'},
    #         $rotamer_angle->{'label_alt_id'};
    #     $rotamer_to_angles{$unique_residue_key}{$rotamer_id}{$rotamer_angle_id}=
    #         $rotamer_angle;
    # }

    # # Determining interaction grid.
    # my $edge_length_interaction =
    #     $parameters->{'_[local]_constants'}{'edge_length_interaction'};
    # my $cutoff_atom =
    #     $parameters->{'_[local]_force_field'}{'cutoff_atom'};

    # # Chooses only CA atoms, because from them, boundary interaction conditions
    # # are measured.
    # my $atom_site_cas =
    #     filter_new( $atom_site,
    #                 { 'include' =>
    #                   { 'label_atom_id' => [ 'CA' ] } } );
    # my ( $grid_box_cas, $grid_ca_atom_pos ) =
    #     grid_box( $parameters, $atom_site_cas, $edge_length_interaction,
    #               extract( $atom_site_cas,
    #                        { 'data' => [ 'id' ], 'is_list' => 1 } ) );
    # my $neighbouring_cells_cas =
    #     identify_neighbour_cells( $grid_box_cas, $grid_ca_atom_pos );

    # for my $grid_id ( keys %{ $grid_ca_atom_pos } ) {
    #     for my $atom_id ( @{ $grid_ca_atom_pos->{$grid_id} } ) {
    #         my $unique_residue_key =
    #             unique_residue_key( $atom_site->{$atom_id} );

    #         $self->{'graph'}->add_vertex( $unique_residue_key );
    #         $self->{'graph'}->set_vertex_attribute(
    #             $unique_residue_key,
    #             'type',
    #             'residue'
    #         );

    #         my $residue_site =
    #             filter_by_unique_residue_key( $atom_site,
    #                                           $unique_residue_key,
    #                                           1 );
    #         connect_atoms( $parameters, $residue_site );
    #         hybridization( $parameters, $residue_site );
    #         rotation_only( $parameters, $residue_site );

    #         $self->{'graph'}->set_vertex_attribute(
    #             $unique_residue_key,
    #             'atom_site',
    #             $residue_site
    #         );

    #         my @rotamer_ids = keys %{ $rotamer_to_angles{$unique_residue_key} };

    #         for my $neighbour_atom_id (@{$neighbouring_cells_cas->{$grid_id}}){
    #             my $neighbour_unique_residue_key =
    #                 unique_residue_key( $atom_site->{$neighbour_atom_id} );

    #             next if $unique_residue_key eq $neighbour_unique_residue_key;

    #             $self->{'graph'}->add_vertex( $neighbour_unique_residue_key );
    #             $self->{'graph'}->set_vertex_attribute(
    #                 $neighbour_unique_residue_key,
    #                 'type',
    #                 'residue'
    #             );

    #             my $neighbour_residue_site =
    #                 filter_by_unique_residue_key( $atom_site,
    #                                               $neighbour_unique_residue_key,
    #                                               1 );
    #             connect_atoms( $parameters, $neighbour_residue_site );
    #             hybridization( $parameters, $neighbour_residue_site );
    #             rotation_only( $parameters, $neighbour_residue_site );

    #             $self->{'graph'}->set_vertex_attribute(
    #                 $neighbour_unique_residue_key,
    #                 'atom_site',
    #                 $neighbour_residue_site
    #             );

    #             next if $self->{'graph'}->has_edge(
    #                 $unique_residue_key, $neighbour_unique_residue_key
    #             );

    #             my $bond_length = bond_length(
    #                 [ [ map { $atom_site->{$atom_id}{$_} }
    #                         ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ],
    #                   [ map { $atom_site->{$neighbour_atom_id}{$_} }
    #                         ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ] ]
    #             );

    #             next if $bond_length > $edge_length_interaction;

    #             $self->{'graph'}->add_edge( $unique_residue_key,
    #                                         $neighbour_unique_residue_key );

    #             my @neighbour_rotamer_ids =
    #                 keys %{ $rotamer_to_angles{$neighbour_unique_residue_key} };

    #             for my $rotamer_id ( @rotamer_ids ) {
    #                 $self->{'graph'}->add_vertex( $rotamer_id );
    #                 $self->{'graph'}->set_vertex_attribute(
    #                     $rotamer_id,
    #                     'type',
    #                     'rotamer'
    #                 );
    #                 $self->{'graph'}->set_vertex_attribute(
    #                     $rotamer_id,
    #                     'angles',
    #                     $rotamer_to_angles{$unique_residue_key}{$rotamer_id}
    #                 );
    #                 $self->{'graph'}->set_vertex_attribute(
    #                     $rotamer_id,
    #                     'energy',
    #                     $rotamer_energies->{$rotamer_id}
    #                 );
    #                 $self->{'graph'}->add_edge( $rotamer_id,
    #                                             $unique_residue_key );


    #                 for my $neighbour_rotamer_id ( @neighbour_rotamer_ids ) {
    #                     $self->{'graph'}->add_vertex( $neighbour_rotamer_id );
    #                     $self->{'graph'}->set_vertex_attribute(
    #                         $neighbour_rotamer_id,
    #                         'type',
    #                         'rotamer'
    #                     );
    #                     $self->{'graph'}->set_vertex_attribute(
    #                         $neighbour_rotamer_id,
    #                         'angles',
    #                         $rotamer_to_angles{$neighbour_unique_residue_key}
    #                                           {$neighbour_rotamer_id}
    #                     );
    #                     $self->{'graph'}->set_vertex_attribute(
    #                         $neighbour_rotamer_id,
    #                         'energy',
    #                         $rotamer_energies->{$neighbour_rotamer_id}
    #                     );
    #                     $self->{'graph'}->add_edge( $neighbour_rotamer_id,
    #                                                 $neighbour_unique_residue_key );
    #                     $self->{'graph'}->add_edge( $rotamer_id,
    #                                                 $neighbour_rotamer_id );
    #                 }
    #             }
    #         }
    #     }
    # }

    return bless $self, $class;
}

# --------------------------------- Methods ---------------------------------- #

sub predict
{
    # my ( $self ) = @_;
    # my $interaction_graph = $self->{'graph'};
}

sub to_tsv
{
    # my ( $self ) = @_;
    # printf "%s\t%s\n", 'node1', 'node2';
    # for my $vertex ( $self->{'graph'}->vertices ) {
    #     for my $neighbour ( $self->{'graph'}->neighbours( $vertex ) ) {
    #         printf "%s\t%s\n", $vertex, $neighbour;
    #     }
    # }
}

# sub predict_sidechains
# {
#     my ( $args ) = @_;
#
#     my %residue_to_atom_site = ();
#     my %rotamer_to_atom_site = ();
#     my %rotamer_pair_energy = ();
#     for my $ca_atom_id ( sort keys %residue_pairs ) {
#         my $unique_residue_key =
#             unique_residue_key( $atom_site->{$ca_atom_id} );

#         if( ! exists $residue_to_atom_site{$unique_residue_key} ) {
#             my $residue_site =
#                 filter_by_unique_residue_key( $atom_site,
#                                               $unique_residue_key,
#                                               1 );
#             connect_atoms( $parameters, $residue_site );
#             hybridization( $parameters, $residue_site );
#             rotation_only( $parameters, $residue_site );
#             $residue_to_atom_site{$unique_residue_key} = $residue_site;
#         }

#         for my $rotamer_id ( keys %{ $residue_to_rotamer{$unique_residue_key}}){
#             my %angles =
#                 map { $rotamer_angles->{$_}{'type'} =>
#                       $rotamer_angles->{$_}{'value'} }
#                     keys %{ $rotamer_to_angles{$rotamer_id} };

#             if( ! exists $rotamer_to_atom_site{$rotamer_id} ) {
#                 my %rotamer_site =
#                     %{ clone( $residue_to_atom_site{$unique_residue_key} ) };
#                 replace_with_rotamer( $parameters, \%rotamer_site,
#                                       $unique_residue_key, \%angles );
#                 $rotamer_to_atom_site{$rotamer_id} = { %rotamer_site };
#             }

#             for my $neighbour_ca_atom_id (
#                 sort keys %{ $residue_pairs{$ca_atom_id} } ) {
#                 my $neighbour_unique_residue_key =
#                     unique_residue_key( $atom_site->{$neighbour_ca_atom_id} );

#                 if( ! exists $residue_to_atom_site{$neighbour_unique_residue_key} ) {
#                     my $neighbour_rotamer_site = filter_by_unique_residue_key(
#                         $atom_site,
#                         $neighbour_unique_residue_key,
#                         1
#                     );
#                     connect_atoms( $parameters, $neighbour_rotamer_site );
#                     hybridization( $parameters, $neighbour_rotamer_site );
#                     rotation_only( $parameters, $neighbour_rotamer_site );
#                     $residue_to_atom_site{$neighbour_unique_residue_key} =
#                         $neighbour_rotamer_site;
#                 }

#                 for my $neighbour_rotamer_id (
#                     keys %{$residue_to_rotamer{$neighbour_unique_residue_key}}){
#                     next if exists $rotamer_pair_energy{$rotamer_id}
#                                                        {$neighbour_rotamer_id} ||
#                             exists $rotamer_pair_energy{$neighbour_rotamer_id}
#                                                        {$rotamer_id};

#                     my @neighbour_angle_ids =
#                         keys %{ $rotamer_to_angles{$neighbour_rotamer_id} };
#                     my %neighbour_angles =
#                         map { $rotamer_angles->{$_}{'type'} =>
#                               $rotamer_angles->{$_}{'value'} }
#                             @neighbour_angle_ids;

#                     if( ! exists $rotamer_to_atom_site{$neighbour_rotamer_id} ){
#                         my %neighbour_rotamer_site =
#                             %{ clone( $residue_to_atom_site{$neighbour_unique_residue_key} ) };
#                         replace_with_rotamer( $parameters,
#                                               \%neighbour_rotamer_site,
#                                               $neighbour_unique_residue_key,
#                                               \%angles );
#                         $rotamer_to_atom_site{$neighbour_rotamer_id} =
#                             { %neighbour_rotamer_site };
#                     }

#                     # Calculate pairwise energy.
#                     my $pairwise_energy_sum = pairwise_rotamer_energy(
#                         $parameters,
#                         $rotamer_to_atom_site{$rotamer_id},
#                         $rotamer_to_atom_site{$neighbour_rotamer_id},
#                         \&ForceField::Bonded::general,
#                         \&ForceField::NonBonded::general,
#                     );

#                     # Under the cutoff limit.
#                     if( $pairwise_energy_sum <= $cutoff_atom  ) {
#                         $rotamer_pair_energy{$rotamer_id}{$neighbour_rotamer_id}=
#                             $pairwise_energy_sum;
#                         $rotamer_pair_energy{$neighbour_rotamer_id}{$rotamer_id}=
#                             $pairwise_energy_sum;
#                     }
#                 }
#             }
#         }
#     }

#     return \%rotamer_pair_energy;
# }

1;
