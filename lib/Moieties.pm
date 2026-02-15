package Moieties;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( missing_atom_names
                     replace_with_moiety
                     replace_hetatoms_with_moiety );
}

use Clone qw( clone );
use List::Util qw( max );
use Math::Trig qw( acos );

use AlterMolecule qw( bond_torsion );
use BondParameters qw( rotatable_bonds );
use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use ForceField::Parameters;
use PDBxParser qw( filter
                   filter_by_unique_residue_key );
use LinearAlgebra qw( mult_matrix_product
                      switch_ref_frame );
use Measure qw( bond_angle );
use PseudoAtoms qw( replace_with_rotamer );
use SidechainModels qw( rotation_translation );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Moieties ---------------------------------- #

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

    $last_atom_id++; # The counter is incresed, because the next will be used.

    my $sig_figs_min = $parameters->{'_[local]_constants'}{'sig_figs_min'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $interaction_atom_names = $parameters->{'_[local]_interaction_atom_names'};
    my $moieties = $parameters->{'_[local]_moieties'};
    my %sidechain_atom_names =
        map { $_ => 1 } @{ $parameters->{'_[local]_sidechain_atom_names'} };
    my $sidechain_hetero_residues =
        $parameters->{'_[local]_sidechain_hetatom_extension'};

    my %all_sidechains = ( %{ $moieties }, %{ $append_moieties } );

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
    my $ca_atom_id =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'CA' ] },
                  'data' => [ 'id' ] } )->[0][0];
    my $ca_atom_coord =
        [ map { $residue_site->{$ca_atom_id}{$_} }
              ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ];
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
        my $moiety_atom = clone $all_sidechains{$moiety}{$atom_id};
        my $moiety_atom_name = $moiety_atom->{'label_atom_id'};
        my $moiety_residue_name = $moiety_atom->{'label_comp_id'};

        next if ! exists $sidechain_atom_names{$moiety_atom_name};

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

        $moiety_atom->{'group_PDB'} =
            exists $sidechain_hetero_residues->{$moiety_residue_name} ?
            'HETATM' : 'ATOM';
        $moiety_atom->{'id'} = $last_atom_id;
        $moiety_atom->{'label_seq_id'} = $residue_id;
        $moiety_atom->{'label_asym_id'} = $residue_chain;
        $moiety_atom->{'auth_seq_id'} =
            $residue_site->{$ca_atom_id}{'auth_seq_id'};
        $moiety_atom->{'auth_asym_id'} =
            $residue_site->{$ca_atom_id}{'auth_asym_id'};
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
        rotatable_bonds( $parameters, $residue_site );

        rotation_translation( $parameters, $residue_site );

        replace_with_rotamer( $parameters, $residue_site, $unique_residue_key,
                              $angles );

        for my $atom_id ( keys %{ $residue_site } ) {
            $atom_site->{$atom_id} = $residue_site->{$atom_id};
        }
    }

    return;
}


#
# Replaces selected heteroatoms with specified moiety (usually another
# heteroatom).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $unique_residue_key - key that can determine unique residue
#     (see PDBxParser::unique_residue_key);
#     $moiety - moiety in atom site data structure.
# Output:
#     changes atom site by replacing heteroatom with specified moiety.
#

sub replace_hetatoms_with_moiety
{
    my ( $parameters, $atom_site, $unique_residue_key, $moiety, $options ) = @_;

    my ( $angles, $append_moieties, $last_atom_id ) =
        ( $options->{'angles'},
          $options->{'append_moieties'},
          $options->{'last_atom_id'}, );

    $angles //= {};
    $append_moieties //= {};
    $last_atom_id //= max( keys %{ $atom_site } );

    $last_atom_id++; # The counter is incresed, because the next will be used.

    my $sig_figs_min = $parameters->{'_[local]_constants'}{'sig_figs_min'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $interaction_atom_names = $parameters->{'_[local]_interaction_atom_names'};
    my $moieties = $parameters->{'_[local]_moieties'};
    my %sidechain_atom_names =
        map { $_ => 1 } @{ $parameters->{'_[local]_sidechain_atom_names'} };
    my $sidechain_hetero_residues =
        $parameters->{'_[local]_sidechain_hetatom_extension'};

    my %all_sidechains = ( %{ $moieties }, %{ $append_moieties } );

    # First, transformation matrix is generated that will position moiety atoms
    # to the origin of reference frame.
    my $moiety_atom_coord =
        filter( { 'atom_site' => $all_sidechains{$moiety},
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my @moiety_helper_z_atom_coord = (
        $moiety_atom_coord->[0],
        $moiety_atom_coord->[1],
        $moiety_atom_coord->[2] + 1,
    );
    my @moiety_helper_y_atom_coord = (
        $moiety_atom_coord->[0],
        $moiety_atom_coord->[1] + 1,
        $moiety_atom_coord->[2],
    );

    my ( $moiety_transf_matrix ) =
        @{ switch_ref_frame( $parameters,
                             $moiety_atom_coord,
                             \@moiety_helper_z_atom_coord,
                             \@moiety_helper_y_atom_coord,
                             'local' ) };

    # Then generates transformation matrix that will align moiety atoms with
    # target atoms.
    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $unique_residue_key, 1 );
    my ( $residue_atom_id ) = sort { $a <=> $b } keys %{ $residue_site };
    my ( $residue_id, $residue_chain, $pdbx_model, $residue_alt,
         $residue_auth_id, $residue_auth_chain, $residue_auth_alt ) =
        split /,/smx, $unique_residue_key;

    my @sidechain_ids =
        @{ filter( { 'atom_site' => $atom_site,
                     'include' =>
                         { 'label_seq_id' => [ $residue_id ],
                           'pdbx_PDB_model_num' => [ $pdbx_model ],
                           'label_asym_id' => [ $residue_chain ] },
                     'data' => [ 'id' ],
                     'is_list' => 1 } ) };

    my $atom_coord =
        filter( { 'atom_site' => $residue_site,
                  'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } )->[0];
    my @helper_z_atom_coord = (
        $atom_coord->[0],
        $atom_coord->[1],
        $atom_coord->[2] + 1,
    );
    my @helper_y_atom_coord = (
        $atom_coord->[0],
        $atom_coord->[1] + 1,
        $atom_coord->[2],
    );

    my ( $transf_matrix ) =
        @{ switch_ref_frame( $parameters,
                             $atom_coord,
                             \@helper_z_atom_coord,
                             \@helper_y_atom_coord,
                             'global' ) };

    # Adds moiety.
    for my $atom_id ( sort keys %{ $all_sidechains{$moiety} } ) {
        my $moiety_atom = clone $all_sidechains{$moiety}{$atom_id};
        my $moiety_atom_name = $moiety_atom->{'label_atom_id'};
        my $moiety_residue_name = $moiety_atom->{'label_comp_id'};

        next if ! exists $sidechain_atom_names{$moiety_atom_name} &&
            ! exists $sidechain_hetero_residues->{$moiety_residue_name};

        my ( $transf_atom_coord ) =
            @{ mult_matrix_product( [ $transf_matrix,
                                      $moiety_transf_matrix,
                                      [ [ $moiety_atom->{'Cartn_x'} ],
                                        [ $moiety_atom->{'Cartn_y'} ],
                                        [ $moiety_atom->{'Cartn_z'} ],
                                        [ 1 ] ] ] ) };

        $moiety_atom->{'group_PDB'} =
            exists $sidechain_hetero_residues->{$moiety_residue_name} ?
            'HETATM' : 'ATOM';
        $moiety_atom->{'id'} = $last_atom_id;
        $moiety_atom->{'label_seq_id'} = $residue_id;
        $moiety_atom->{'label_asym_id'} = $residue_chain;
        $moiety_atom->{'label_entity_id'} =
            $residue_site->{$residue_atom_id}{'label_entity_id'};
        $moiety_atom->{'auth_seq_id'} =
            $residue_site->{$residue_atom_id}{'auth_seq_id'};
        $moiety_atom->{'auth_asym_id'} =
            $residue_site->{$residue_atom_id}{'auth_asym_id'};
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
        rotatable_bonds( $parameters, $residue_site );

        rotation_translation( $parameters, $residue_site );

        replace_with_rotamer( $parameters, $residue_site, $unique_residue_key,
                              $angles );

        for my $atom_id ( keys %{ $residue_site } ) {
            $atom_site->{$atom_id} = $residue_site->{$atom_id};
        }
    }

    return;
}

#
# Checks if there are all mandatory atoms in the residue.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $unique_residue_key - key that can determine unique residue
#     (see PDBxParser::unique_residue_key).
# Output:
#     outputs true or false.
#

sub missing_atom_names
{
    my ( $parameters, $residue_site ) = @_;

    my @atom_ids = keys %{ $residue_site };
    my $residue_name = $residue_site->{$atom_ids[0]}{'label_comp_id'};

    my %mandatory_residue_atoms =
        map { $_ => 0 }
        keys %{ $parameters->{'_[local]_residue_atom_necessity'}{$residue_name}{'mandatory'} };

    for my $atom_name ( map { $residue_site->{$_}{'label_atom_id'} }
                        keys %{ $residue_site } ) {
        $mandatory_residue_atoms{$atom_name} = 1;
    }

    my @missing_atom_names =
        sort
        grep { ! $mandatory_residue_atoms{$_} }
        keys %mandatory_residue_atoms;

    return \@missing_atom_names;
}

1;
