package Moieties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %atoms
                     %sidechain
                     replace_with_moiety );

use List::Util qw( max );
use Math::Trig qw( acos );

use AlterMolecule qw( bond_torsion );
use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use ForceField::Parameters;
use PDBxParser qw( filter
                   filter_by_unique_residue_key );
use LinearAlgebra qw( mult_matrix_product
                      switch_ref_frame );
use Measure qw( bond_angle );
use PseudoAtoms qw( replace_with_rotamer );
use SidechainModels qw( rotation_only );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Moieties ---------------------------------- #

our %atoms = (
    'H' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'H',
            'label_alt_id' => q{.},
            'Cartn_x' => -53.464,
            'Cartn_y' => -109.890,
            'Cartn_z' => -5.544,
            'pdbx_PDB_model_num' => 1,
        },
    },
);

our %sidechains = (
    'SER' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.009,
            'Cartn_y' => 0.077,
            'Cartn_z' => -0.688,
            'pdbx_PDB_model_num' => 1,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.494,
            'Cartn_y' => 0.929,
            'Cartn_z' => 0.504,
            'pdbx_PDB_model_num' => 1,
        },
        3 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'O',
            'label_atom_id' => 'OG',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.029,
            'Cartn_y' => 0.446,
            'Cartn_z' => 1.769,
            'pdbx_PDB_model_num' => 1,
        },
    },
    'PHE' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 1.893,
            'Cartn_y' => -0.429,
            'Cartn_z' => -1.043,
            'pdbx_PDB_model_num' => 1,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 1.385,
            'Cartn_y' => 0.334,
            'Cartn_z' => 0.218,
            'pdbx_PDB_model_num' => 1,
        },
        3 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'C',
            'label_atom_id' => 'CG',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.139,
            'Cartn_y' => 0.429,
            'Cartn_z' => 0.393,
            'pdbx_PDB_model_num' => 1,
        },
        4 => {
            'group_PDB' => 'ATOM',
            'id' => 4,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.821,
            'Cartn_y' => -0.586,
            'Cartn_z' => 1.074,
            'pdbx_PDB_model_num' => 1,
        },
        5 => {
            'group_PDB' => 'ATOM',
            'id' => 5,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.860,
            'Cartn_y' => 1.483,
            'Cartn_z' => -0.173,
            'pdbx_PDB_model_num' => 1,
        },
        6 => {
            'group_PDB' => 'ATOM',
            'id' => 6,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.208,
            'Cartn_y' => -0.556,
            'Cartn_z' => 1.171,
            'pdbx_PDB_model_num' => 1,
        },
        7 => {
            'group_PDB' => 'ATOM',
            'id' => 7,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.248,
            'Cartn_y' => 1.518,
            'Cartn_z' => -0.066,
            'pdbx_PDB_model_num' => 1,
        },
        8 => {
            'group_PDB' => 'ATOM',
            'id' => 8,
            'type_symbol' => 'C',
            'label_atom_id' => 'CZ',
            'label_alt_id' => q{.},
            'label_comp_id' => 'PHE',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.921,
            'Cartn_y' => 0.497,
            'Cartn_z' => 0.603,
            'pdbx_PDB_model_num' => 1,
        },
    },
    'TRP' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => ,
            'Cartn_y' => ,
            'Cartn_z' => ,
            'pdbx_PDB_model_num' => 1,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => ,
            'Cartn_y' => ,
            'Cartn_z' => ,
            'pdbx_PDB_model_num' => 1,
        },
        3 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'C',
            'label_atom_id' => 'CG',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => ,
            'Cartn_y' => ,
            'Cartn_z' => ,
            'pdbx_PDB_model_num' => 1,
        },
        4 => {
            'group_PDB' => 'ATOM',
            'id' => 4,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => ,
            'Cartn_y' => ,
            'Cartn_z' => ,
            'pdbx_PDB_model_num' => 1,
        },
        5 => {
            'group_PDB' => 'ATOM',
            'id' => 5,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.515,
            'Cartn_y' => 0.672,
            'Cartn_z' => 0.130,
            'pdbx_PDB_model_num' => 1,
        },
        6 => {
            'group_PDB' => 'ATOM',
            'id' => 6,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -1.895,
            'Cartn_y' => 0.486,
            'Cartn_z' => 0.393,
            'pdbx_PDB_model_num' => 1,
        },
        7 => {
            'group_PDB' => 'ATOM',
            'id' => 7,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE3',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.079,
            'Cartn_y' => 1.700,
            'Cartn_z' => -0.745,
            'pdbx_PDB_model_num' => 1,
        },
        8 => {
            'group_PDB' => 'ATOM',
            'id' => 8,
            'type_symbol' => 'N',
            'label_atom_id' => 'NE1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.092,
            'Cartn_y' => -0.588,
            'Cartn_z' => 1.248,
            'pdbx_PDB_model_num' => 1,
        },
        9 => {
            'group_PDB' => 'ATOM',
            'id' => 9,
            'type_symbol' => 'C',
            'label_atom_id' => 'CZ2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.849,
            'Cartn_y' => 1.335,
            'Cartn_z' => -0.209,
            'pdbx_PDB_model_num' => 1,
        },
        10 => {
            'group_PDB' => 'ATOM',
            'id' => 10,
            'type_symbol' => 'C',
            'label_atom_id' => 'CZ3',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -1.040,
            'Cartn_y' => 2.530,
            'Cartn_z' => -1.322,
            'pdbx_PDB_model_num' => 1,
        },
        11 => {
            'group_PDB' => 'ATOM',
            'id' => 11,
            'type_symbol' => 'C',
            'label_atom_id' => 'CH2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TRP',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -2.405,
            'Cartn_y' => 2.353,
            'Cartn_z' => -1.057,
            'pdbx_PDB_model_num' => 1,
        },
    },
    'TYR' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.612,
            'Cartn_y' => -2.082,
            'Cartn_z' => 1.127,
            'pdbx_PDB_model_num' => 1,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 1.052,
            'Cartn_y' => -0.589,
            'Cartn_z' => 1.096,
            'pdbx_PDB_model_num' => 1,
        },
        3 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'C',
            'label_atom_id' => 'CG',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.390,
            'Cartn_y' => 0.302,
            'Cartn_z' => 0.038,
            'pdbx_PDB_model_num' => 1,
        },
        4 => {
            'group_PDB' => 'ATOM',
            'id' => 4,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.807,
            'Cartn_y' => 0.254,
            'Cartn_z' => -1.296,
            'pdbx_PDB_model_num' => 1,
        },
        5 => {
            'group_PDB' => 'ATOM',
            'id' => 5,
            'type_symbol' => 'C',
            'label_atom_id' => 'CD2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.647,
            'Cartn_y' => 1.164,
            'Cartn_z' => 0.406,
            'pdbx_PDB_model_num' => 1,
        },
        6 => {
            'group_PDB' => 'ATOM',
            'id' => 6,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE1',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.188,
            'Cartn_y' => 1.057,
            'Cartn_z' => -2.250,
            'pdbx_PDB_model_num' => 1,
        },
        7 => {
            'group_PDB' => 'ATOM',
            'id' => 7,
            'type_symbol' => 'C',
            'label_atom_id' => 'CE2',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -1.264,
            'Cartn_y' => 1.966,
            'Cartn_z' => -0.549,
            'pdbx_PDB_model_num' => 1,
        },
        8 => {
            'group_PDB' => 'ATOM',
            'id' => 8,
            'type_symbol' => 'C',
            'label_atom_id' => 'CZ',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.846,
            'Cartn_y' => 1.912,
            'Cartn_z' => -1.877,
            'pdbx_PDB_model_num' => 1,
        },
        9 => {
            'group_PDB' => 'ATOM',
            'id' => 9,
            'type_symbol' => 'O',
            'label_atom_id' => 'OH',
            'label_alt_id' => q{.},
            'label_comp_id' => 'TYR',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -1.452,
            'Cartn_y' => 2.696,
            'Cartn_z' => -2.817,
            'pdbx_PDB_model_num' => 1,
        },
    },
);

#
# Replaces selected side-chain with specified moiety (usually another
# side-chain).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $unique_residue_key - key that can determine unique residue
#     (see PDBxParser::unique_residue_key);
#     $moiety - moiety in atom site data structure.
# Output:
#     changes atom site by replacing side-chain with specified moiety.
#

sub replace_with_moiety
{
    my ( $parameters, $atom_site, $unique_residue_key, $moiety, $options ) = @_;

    my ( $isomer, $angles, $append_moieties, $last_atom_id ) =
        ( $options->{'isomer'},
          $options->{'angles'},
          $options->{'append_moieties'},
          $options->{'last_atom_id'}, );

    $isomer //= 'R';
    $angles //= {};
    $append_moieties //= {};
    $last_atom_id //= max( keys %{ $atom_site } );

    my $sig_figs_min = $parameters->{'_[local]_constants'}{'sig_figs_min'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $interaction_atom_names = $parameters->{'_[local]_interaction_atom_names'};

    my %all_sidechains = ( %sidechains, %{ $append_moieties } );

    # First, transformation matrix is generated that will position moiety atoms
    # to the origin of reference frame.
    my $moiety_ca_atom_coord =
        filter( { 'atom_site' => $all_sidechains{$moiety},
                  'include' => { 'label_atom_id' => [ 'CA' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my $moiety_cb_atom_coord =
        filter( { 'atom_site' => $all_sidechains{$moiety},
                  'include' => { 'label_atom_id' => [ 'CB' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my @moiety_helper_atom_coord = map { $_ + 1 } @{ $moiety_ca_atom_coord };

    my ( $moiety_transf_matrix ) =
        @{ switch_ref_frame( $parameters,
                             $moiety_ca_atom_coord,
                             $moiety_cb_atom_coord,
                             \@moiety_helper_atom_coord,
                             'local' ) };

    # Then generates transformation matrix that will align moiety atoms with
    # target atoms.
    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $unique_residue_key, 1 );
    my ( $residue_id, $residue_chain, $pdbx_model, $residue_alt ) =
        split /,/smx, $unique_residue_key;

    my @sidechain_ids =
        @{ filter( { 'atom_site' => $atom_site,
                     'include' =>
                         { 'label_seq_id' => [ $residue_id ],
                           'pdbx_PDB_model_num' => [ $pdbx_model ],
                           'label_asym_id' => [ $residue_chain ] },
                     'exclude' =>
                         # TODO: make proper list of mainchain atoms.
                         { 'label_atom_id' =>
                               [ grep { $_ ne 'CB' }
                                      @{ $interaction_atom_names } ] },
                     'data' => [ 'id' ],
                     'is_list' => 1 } ) };

    my $n_atom_coord =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'N' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my $ca_atom_coord =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'CA' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my $c_atom_coord =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'C' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my $o_atom_coord =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'O' ] },
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];

    # TODO: should be refactored, because the code is familiar to
    # PseudoAtoms::add_hydrogen().
    my $bond_angle =
        bond_angle( [ $n_atom_coord, $ca_atom_coord, $c_atom_coord ] );

    my $moiety_angle = acos( ( - 4 - 2 * cos $bond_angle ) / 10 );

    my ( $transf_matrix ) =
        @{ switch_ref_frame( $parameters,
                             $ca_atom_coord,
                             $n_atom_coord,
                             $c_atom_coord,
                             'global' ) };

    # Rotational matrix is created for producing 'R' or 'S' configuration.
    my $rotational_matrix = bond_torsion( $parameters,
                                          $c_atom_coord,
                                          $ca_atom_coord,
                                          $o_atom_coord,
                                          'omega' );

    # Adds moiety.
    for my $atom_id ( sort keys %{ $all_sidechains{$moiety} } ) {
        my $moiety_atom = $all_sidechains{$moiety}{$atom_id};
        my ( $transf_atom_coord ) =
            @{ mult_matrix_product( [ @{ $rotational_matrix },
                                      $transf_matrix,
                                      $moiety_transf_matrix,
                                      [ [ $moiety_atom->{'Cartn_x'} ],
                                        [ $moiety_atom->{'Cartn_y'} ],
                                        [ $moiety_atom->{'Cartn_z'} ],
                                        [ 1 ] ] ],
                                    { $isomer eq 'R' ?
                                      ( 'omega' => (-2) * $pi / 3 ) :
                                      ( 'omega' =>   2  * $pi / 3 ) } ) };

        $moiety_atom->{'id'} = $last_atom_id;
        $moiety_atom->{'label_seq_id'} = $residue_id;
        $moiety_atom->{'label_asym_id'} = $residue_chain;
        $moiety_atom->{'pdbx_PDB_model_num'} = $pdbx_model;
        # TODO: check if there will be situations when non '.' label_alt_id is
        # needed.
        $moiety_atom->{'label_alt_id'} = q{.};
        $moiety_atom->{'Cartn_x'} =
            sprintf $sig_figs_min, $transf_atom_coord->[0][0];
        $moiety_atom->{'Cartn_y'}=
            sprintf $sig_figs_min, $transf_atom_coord->[1][0];
        $moiety_atom->{'Cartn_z'}=
            sprintf $sig_figs_min, $transf_atom_coord->[2][0];

        $atom_site->{$last_atom_id} = $moiety_atom;
        $last_atom_id++;
    }

    # Removes old side-chain atoms.
    foreach( @sidechain_ids ) {
        delete $residue_site->{$_};
        delete $atom_site->{$_};
    }

    # Renames residue.
    foreach( keys %{ $residue_site } ) {
        $atom_site->{$_}{'label_comp_id'} = $moiety;
    }

    if( %{ $angles } ) {
        $residue_site =
            filter_by_unique_residue_key( $atom_site, $unique_residue_key, 1 );

        connect_atoms( $parameters, $residue_site );
        hybridization( $parameters, $residue_site );

        rotation_only( $parameters, $residue_site );

        replace_with_rotamer( $parameters, $residue_site, $unique_residue_key,
                              $angles );

        for my $atom_id ( keys %{ $residue_site } ) {
            $atom_site->{$atom_id} = $residue_site->{$atom_id};
        }
    }

    return;
}

1;
