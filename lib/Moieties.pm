package Moieties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %ATOMS
                     %SIDECHAINS
                     replace_with_moiety );

use List::Util qw( max );
use Math::Trig qw( acos );

use AlterMolecule qw( bond_torsion );
use AtomProperties qw( @MAINCHAIN_NAMES );
use Constants qw( $PI
                  $SIG_FIGS );
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

our %ATOMS = (
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

our %SIDECHAINS = (
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
    my ( $atom_site, $unique_residue_key, $moiety, $options ) = @_;

    my ( $isomer, $angles, $append_moieties, $last_atom_id ) =
        ( $options->{'isomer'},
          $options->{'angles'},
          $options->{'append_moieties'},
          $options->{'last_atom_id'}, );

    $isomer //= 'R';
    $angles //= {};
    $append_moieties //= {};
    $last_atom_id //= max( keys %{ $atom_site } );

    my %all_sidechains = ( %SIDECHAINS, %{ $append_moieties } );

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
        @{ switch_ref_frame( $moiety_ca_atom_coord,
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
                               [ grep { $_ ne 'CB' } @MAINCHAIN_NAMES ] },
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
        @{ switch_ref_frame( $ca_atom_coord,
                             $n_atom_coord,
                             $c_atom_coord,
                             'global' ) };

    # Rotational matrix is created for producing 'R' or 'S' configuration.
    my $rotational_matrix = bond_torsion( $c_atom_coord,
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
                                      ( 'omega' => (-2) * $PI / 3 ) :
                                      ( 'omega' =>   2  * $PI / 3 ) } ) };

        $moiety_atom->{'id'} = $last_atom_id;
        $moiety_atom->{'label_seq_id'} = $residue_id;
        $moiety_atom->{'label_asym_id'} = $residue_chain;
        $moiety_atom->{'pdbx_PDB_model_num'} = $pdbx_model;
        # TODO: check if there will be situations when non '.' label_alt_id is
        # needed.
        $moiety_atom->{'label_alt_id'} = q{.};
        $moiety_atom->{'Cartn_x'}= sprintf $SIG_FIGS, $transf_atom_coord->[0][0];
        $moiety_atom->{'Cartn_y'}= sprintf $SIG_FIGS, $transf_atom_coord->[1][0];
        $moiety_atom->{'Cartn_z'}= sprintf $SIG_FIGS, $transf_atom_coord->[2][0];

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

        rotation_only( $residue_site );
        replace_with_rotamer( $residue_site, $unique_residue_key, $angles );

        for my $atom_id ( keys %{ $residue_site } ) {
            $atom_site->{$atom_id} = $residue_site->{$atom_id};
        }
    }

    return;
}

1;
