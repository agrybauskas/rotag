package Hydrogens;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( add_hydrogens );

use List::Util qw( max );
use List::MoreUtils qw( any );
use Math::Trig qw( acos );

use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use LinearAlgebra qw( matrix_product
                      mult_matrix_product
                      switch_ref_frame );
use Measure qw( bond_angle
                dihedral_angle );
use PDBxParser qw( create_pdbx_entry
                   filter );

# -------------------------- Generation of hydrogen atoms --------------------- #

#
# Adds hydrogens to the molecule.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry. For example, hydrogens of
#     methyl group that has only one known connection are not added, but with two
#     connections - are;
#     $options->{use_existing_connections} - does not overwrite atom connections;
#     $options->{use_existing_hybridizations} - does not overwrite atom
#     hybridizations;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections;
#     $options->{exclude_by_atom_name} - list of atom names that are excluded
#     from being added hydrogens to;
#     $options->{alt_group_id} - alternative group id that is used to distinguish
#     pseudo atoms.
# Output:
#     %hydrogen_site - atom data structure with added hydrogens.
#

sub add_hydrogens
{
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $add_only_clear_positions, $use_existing_connections,
         $use_existing_hybridizations, $reference_atom_site,
         $exclude_by_atom_name, $exclude_by_atom_ids,
         $last_atom_id, $alt_group_id,
         $use_origins_alt_group_id, $selection_state ) =
        ( $options->{'add_only_clear_positions'},
          $options->{'use_existing_connections'},
          $options->{'use_existing_hybridizations'},
          $options->{'reference_atom_site'},
          $options->{'exclude_by_atom_name'},
          $options->{'exclude_by_atom_ids'},
          $options->{'last_atom_id'},
          $options->{'alt_group_id'},
          $options->{'use_origins_alt_group_id'},
          $options->{'selection_state'}, );

    $add_only_clear_positions //= 0;
    $use_existing_connections //= 0;
    $use_existing_hybridizations //= 0;
    $reference_atom_site //= $atom_site;
    $exclude_by_atom_name //= [];
    $exclude_by_atom_ids //= [];
    $alt_group_id //= '1';
    $use_origins_alt_group_id //= 0;

    my $hydrogen_names = $parameters->{'_[local]_hydrogen_names'};
    my $residue_atoms = $parameters->{'_[local]_residue_atom_necessity'};
    my $connectivity = $parameters->{'_[local]_connectivity'};
    my $sig_figs_min = $parameters->{'_[local]_constants'}{'sig_figs_min'};

    my %atom_site = %{ $atom_site };

    if( !$use_existing_connections ) {connect_atoms( $parameters, \%atom_site )};

    if( !$use_existing_hybridizations ) {hybridization($parameters,\%atom_site,)};

    my %hydrogen_site;
    $last_atom_id //= max( keys %{ $atom_site } );

    for my $atom_id ( sort { $a <=> $b } keys %atom_site ) {
        my $atom_name = $atom_site{$atom_id}{'label_atom_id'};

        next if any { $_ eq $atom_id } @{ $exclude_by_atom_ids };
        next if any { $_ eq $atom_name } @{ $exclude_by_atom_name };

        my $residue_name = $atom_site{$atom_id}{'label_comp_id'};

        my $hydrogen_names = $hydrogen_names->{$residue_name}{$atom_name};

        if( ! $hydrogen_names ) { next; }; # Exits early if there should be no
                                           # hydrogens connected to the atom.

        my $hybridization = $atom_site->{$atom_id}{'hybridization'};

        # Decides how many and what hydrogens should be added according to the
        # quantity of bonds and hydrogen atoms that should be connected to the
        # target atom.
        my @connection_ids = ();
        if( exists $reference_atom_site->{"$atom_id"}{'connections'} ) {
            @connection_ids =
                @{ $reference_atom_site->{"$atom_id"}{'connections'} };
        }

        # TODO: should be pre-determined as constant variable.
        my @mandatory_residue_atoms =
            keys %{ $residue_atoms->{$residue_name}{'mandatory'} };
        my @mandatory_connections = ();
        for my $mandatory_atom ( @mandatory_residue_atoms ) {
            if( any { $mandatory_atom eq $_ }
                   @{ $connectivity->{$residue_name}{$atom_name} } ) {
                push @mandatory_connections, $mandatory_atom;
            }
        }

        # Hydrogens cannot be added if there is a missing information about
        # mandatory atom connections.
        next if( scalar @connection_ids < scalar @mandatory_connections );

        my @connection_names =
            map { $reference_atom_site->{"$_"}{'label_atom_id'} } @connection_ids;
        my @missing_hydrogens;
        for my $hydrogen_name ( @{ $hydrogen_names } ) {
            if( ! any { /$hydrogen_name/sxm } @connection_names ) {
                push @missing_hydrogens, $hydrogen_name;
            }
        }

        # Exits early if there are no spare hydrogens to add.
        if( scalar @missing_hydrogens == 0 ) { next; };

        #            sp3                       sp2               sp
        #
        #            Up(2)                     Up(2)             Up(2)
        # z          |                         |                 |
        # |_y      Middle(1) __ Right(3)     Middle(1)         Middle(1)
        # /         / \                       / \                |
        # x    Left(4) Back(5)           Left(4) Right(3)        Down(3)
        #
        # Depending on hybridization and present bond connections, adds missing
        # hydrogens.
        my %hydrogen_coord = map { $_ => undef } @missing_hydrogens;

        if( $hybridization eq 'sp3' ) {
            add_hydrogens_sp3( $parameters, $atom_site, $atom_id,
                               \%hydrogen_coord, \@missing_hydrogens, $options );
        } elsif( $hybridization eq 'sp2' ) {
            add_hydrogens_sp2( $parameters, $atom_site, $atom_id,
                               \%hydrogen_coord, \@missing_hydrogens, $options );
        } elsif( $hybridization eq 'sp' ) {
            add_hydrogens_sp( $parameters, $atom_site, $atom_id,
                              \%hydrogen_coord, \@missing_hydrogens, $options );
        }

        # Each coordinate of atoms is transformed by transformation
        # matrix and added to %hydrogen_site.
        for my $hydrogen_name ( sort { $a cmp $b } keys %hydrogen_coord ) {
            if( $hydrogen_coord{$hydrogen_name} ) {
            # Adds necessary PDBx entries to pseudo atom site.
            $last_atom_id++;
            create_pdbx_entry(
                { 'group_PDB' => $atom_site->{$atom_id}{'group_PDB'},
                  'atom_site' => \%hydrogen_site,
                  'id' => $last_atom_id,
                  'type_symbol' => 'H',
                  'label_atom_id' => $hydrogen_name,
                  'label_alt_id' =>
                      $use_origins_alt_group_id ?
                      $atom_site->{$atom_id}{'label_alt_id'} : $alt_group_id,
                  'label_comp_id' => $atom_site->{$atom_id}{'label_comp_id'},
                  'label_asym_id' => $atom_site->{$atom_id}{'label_asym_id'},
                  'label_entity_id' => $atom_site->{$atom_id}{'label_entity_id'},
                  'label_seq_id' => $atom_site{$atom_id}{'label_seq_id'},
                  'cartn_x' =>
                      sprintf( $sig_figs_min,
                               $hydrogen_coord{$hydrogen_name}->[0][0] ),
                  'cartn_y' =>
                      sprintf( $sig_figs_min,
                               $hydrogen_coord{$hydrogen_name}->[1][0] ),
                  'cartn_z' =>
                      sprintf( $sig_figs_min,
                               $hydrogen_coord{$hydrogen_name}->[2][0] ),
                  'pdbx_PDB_model_num' =>
                      $atom_site->{$atom_id}{'pdbx_PDB_model_num'},
                } );
            # Adds additional pseudo-atom flag for future filtering.
            $hydrogen_site{$last_atom_id}{'is_pseudo_atom'} = 1;
            # Adds atom id that pseudo atoms was made of.
            $hydrogen_site{$last_atom_id}{'origin_atom_id'} = $atom_id;
            # Marks origin atom id as connection.
            $hydrogen_site{$last_atom_id}{'connections'} = [ $atom_id ];
            # By default, hydrogen atoms are sp3 hybridized.
            $hydrogen_site{$last_atom_id}{'hybridization'} = 'sp3';
            # Adds selection state if it is defined.
            if( defined $selection_state ) {
                $hydrogen_site{$last_atom_id}{'[local]_selection_state'} =
                    $selection_state;
            } else {
                $hydrogen_site{$last_atom_id}{'[local]_selection_state'} =
                    $atom_site->{$atom_id}{'[local]_selection_state'};
            }
            }
        }
    }

    return \%hydrogen_site;
}

#
# Adds hydrogens to the atoms which are sp3 hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp3
{
    my ( $parameters, $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;

    my ( $add_only_clear_positions, $reference_atom_site ) = (
        $options->{'add_only_clear_positions'},
        $options->{'reference_atom_site'},
    );

    my $atom_properties = $parameters->{'_[local]_atom_properties'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    $add_only_clear_positions //= 0;
    $reference_atom_site //= $atom_site;

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $atom_properties->{$atom_type}{'covalent_radius'}{'length'}[0] +
        $atom_properties->{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids =
        exists $reference_atom_site->{"$atom_id"}{'connections'} ?
        @{ $reference_atom_site->{"$atom_id"}{'connections'} } :
        ();
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my $lone_pair_count = $atom_properties->{$atom_type}{'lone_pairs'};

    if( scalar( @connection_ids ) == 3 ) {
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord,
             $right_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]},
                   $atom_coord{$connection_ids[2]}, );

        # Theoretically optimal angles should be so, the distance of
        # atoms would be furthest from each other. This strategy could
        # be achieved by imagining each bond as vector of the length 1.
        # Then the sum of all vectors would be 0, if angles are equal.
        #
        #                        ->  ->  ->  ->
        #                        A + B + C + D = 0
        #
        # In that case, the square of sum also should be equal 0.
        #
        #                        ->  ->  ->  ->
        #                      ( A + B + C + D ) ^ 2 = 0
        #
        #       ->  ->      ->  ->      ->  ->      ->  ->      ->  ->
        #   2 * A * B + 2 * A * C + 2 * A * D + 2 * B * C + 2 * B * D +
        #       ->  ->   ->      ->      ->      ->
        # + 2 * C + D  + A ^ 2 + B ^ 2 + C ^ 2 + D ^ 2 = 0
        #
        # Because length of each vector is equal to 1, then dot product
        # is equal to cos(alpha), where all angles between bonds are
        # equal. If there were no given angles, alpha should be 109.5
        # degrees.
        #
        #                   alpha = arccos( - 1 / 3 )
        #
        # However, there is a restriction of given angle. And the
        # calculation changes:
        #
        #        alpha = arccos( ( - 4 - 2 * cos( beta ) -
        #                                2 * cos( gamma ) -
        #                                2 * cos( delta ) ) / 6 )
        #
        # where beta is the given angle.

        # Calculates all bond angles between atoms connected to the
        # middle atom.
        my $up_mid_right_angle =
            bond_angle( [ $up_atom_coord,
                          $mid_atom_coord,
                          $right_atom_coord ] );
        my $up_mid_left_angle =
            bond_angle( [ $up_atom_coord,
                          $mid_atom_coord,
                          $left_atom_coord ] );
        my $right_mid_left_angle =
            bond_angle( [ $right_atom_coord,
                          $mid_atom_coord,
                          $left_atom_coord ] );

        # Calculates what common angle between hydrogen and rest of the
        # atoms should be.
        my $hydrogen_angle =
            acos( ( - 4 -
                    2 * cos( $up_mid_right_angle ) -
                    2 * cos( $up_mid_left_angle ) -
                    2 * cos $right_mid_left_angle ) /
                  6 );

        # Determines dihedral angle between left and right atoms. Then
        # splits rest of the 2 * pi angle into two equal parts.
        my $dihedral_angle =
            dihedral_angle( [ $left_atom_coord,
                              $up_atom_coord,
                              $mid_atom_coord,
                              $right_atom_coord ] );
        if( abs( $dihedral_angle ) < ( 3 * $pi / 4 ) ) {
            if( $dihedral_angle < 0 ) {
                $dihedral_angle =   ( 2 * $pi + $dihedral_angle ) / 2;
            } else {
                $dihedral_angle = - ( 2 * $pi - $dihedral_angle ) / 2;
            }
        } else {
            if( $dihedral_angle < 0 ) {
                $dihedral_angle =   $dihedral_angle / 2;
            } else {
                $dihedral_angle = - $dihedral_angle / 2;
            }
        }

        # Places hydrogen according to previously calculated angles.
        my ( $transf_matrix ) =
            @{ switch_ref_frame(
                   $parameters,
                   $mid_atom_coord,
                   $up_atom_coord,
                   $left_atom_coord,
                   'global' ) };

        ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
            @{ mult_matrix_product(
                   [ $transf_matrix,
                     [ [ $bond_length *
                         cos( $pi / 2 - $dihedral_angle ) *
                         sin $hydrogen_angle ],
                       [ $bond_length *
                         sin( $pi / 2 - $dihedral_angle ) *
                         sin $hydrogen_angle ],
                       [ $bond_length *
                         cos $hydrogen_angle ],
                       [ 1 ] ] ] ) };
    } elsif( scalar( @connection_ids ) == 2 &&
           ! ( $add_only_clear_positions && $lone_pair_count > 0 ) ) {
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]}, );

        my $bond_angle =
            bond_angle( [ $up_atom_coord,
                          $mid_atom_coord,
                          $left_atom_coord ] );

        # This time is only one defined angle. And the calculation
        # changes:
        #
        #       acos( ( - 4 - 2 * cos( $bond_angle ) ) / 10 )
        #
        my $hydrogen_angle = acos( ( - 4 - 2 * cos $bond_angle ) / 10 );

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $left_atom_coord,
                                 'global' ) };

        # Adds hydrogen first to both atoms that have 0 or 1 electron
        # pairs.
        if( scalar @{ $missing_hydrogens } >= 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 7 * $pi / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             sin( 7 * $pi / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             cos $hydrogen_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        # Additional hydrogen is added only to the atom that has no
        # electron pairs.
        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( - 1 * $pi / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             sin( - 1 * $pi / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             cos $hydrogen_angle ],
                           [ 1 ] ] ] ) };
        }

    } elsif( scalar @connection_ids == 1 && ! $add_only_clear_positions ) {
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $side_coord ) = # Coordinate that only will be used for
                 # defining a local reference frame.
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   [ $atom_site->{$atom_id}{'Cartn_x'},
                     $atom_site->{$atom_id}{'Cartn_y'} + 1,
                     $atom_site->{$atom_id}{'Cartn_z'}, ], );

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $side_coord,
                                 'global' ) };

        # Decreases bond angle, if lone pairs are present.
        # TODO: check if angle reduction is relevant in amino acid
        # structures.
        my $bond_angle;
        if( $lone_pair_count > 0 ) {
            $bond_angle = ( 109.5 - $lone_pair_count * 2.5 ) * $pi / 180;
        } else {
            $bond_angle = 109.5 * $pi / 180;
        }

        # Adds hydrogens according to the quantity of lone pairs.
        if( scalar @{ $missing_hydrogens } >= 3 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length * sin $bond_angle ],
                           [ 0 ],
                           [ $bond_length * cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } >= 2 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 2 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 2 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 4 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 4 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
        }
    } elsif( scalar @connection_ids == 0 && ! $add_only_clear_positions ) {
        # TODO: should be refactored as similar instances of code can be
        # observed.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $side_coord ) =
            ( [ $atom_site->{$atom_id}{'Cartn_x'},
                $atom_site->{$atom_id}{'Cartn_y'},
                $atom_site->{$atom_id}{'Cartn_z'} + 1, ],
              $atom_coord{$atom_id},
              [ $atom_site->{$atom_id}{'Cartn_x'},
                $atom_site->{$atom_id}{'Cartn_y'} + 1,
                $atom_site->{$atom_id}{'Cartn_z'}, ], );

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $side_coord,
                                 'global' ) };

        my $bond_angle;
        if( $lone_pair_count > 0 ) {
            $bond_angle = ( 109.5 - $lone_pair_count * 2.5 ) * $pi / 180;
        } else {
            $bond_angle = 109.5 * $pi / 180;
        }

        if( scalar @{ $missing_hydrogens } >= 3 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length * sin $bond_angle ],
                           [ 0 ],
                           [ $bond_length * cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } >= 2 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 2 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 2 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 4 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 4 * $pi / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
        }
    }

    return;
}

#
# Adds hydrogens to the atoms which are sp2 hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp2
{
    my ( $parameters, $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;

    my ( $add_only_clear_positions, $reference_atom_site ) = (
        $options->{'add_only_clear_positions'},
        $options->{'reference_atom_site'}
    );

    my $atom_properties = $parameters->{'_[local]_atom_properties'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    $add_only_clear_positions //= 0;
    $reference_atom_site //= $atom_site;

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $atom_properties->{$atom_type}{'covalent_radius'}{'length'}[1] +
        $atom_properties->{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids = @{ $reference_atom_site->{"$atom_id"}{'connections'} };
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my $lone_pair_count = $atom_properties->{$atom_type}{'lone_pairs'};

    # Depending on quantity of atoms connections, adds hydrogens.
    if( scalar @connection_ids == 2 ) {
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]}, );

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $left_atom_coord,
                                 'global' ) };

        # Calculates angle between two bonds where target atom takes
        # part. Then desirable hydrogen angle should be calculated by
        # formula:
        #
        #                  angle = ( 360 - beta ) / 2
        #
        # where beta is angle between two given bonds.
        my $bond_angle =
            ( 2 * $pi - bond_angle( [ $up_atom_coord,
                                      $mid_atom_coord,
                                      $left_atom_coord ] ) ) / 2;

        # Hydrogen is placed by placing hydrogen colinearly and then
        # rotating according to bond angle.
        ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
            @{ mult_matrix_product(
                   [ $transf_matrix,
                     [ [ $bond_length *
                         cos( - 0.5 * $pi ) *
                         sin $bond_angle ],
                       [ $bond_length *
                         sin( - 0.5 * $pi ) *
                         sin $bond_angle ],
                       [ $bond_length *
                         cos $bond_angle ],
                       [ 1 ] ] ] ) };

    } elsif( ( scalar @connection_ids == 1 ) &&
           ! ( $add_only_clear_positions && $lone_pair_count > 1 ) ) {
        my ( $up_atom_coord,
             $mid_atom_coord,
             $side_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   [ $atom_site->{$atom_id}{'Cartn_x'},
                     $atom_site->{$atom_id}{'Cartn_y'} + 1,
                     $atom_site->{$atom_id}{'Cartn_z'}, ], );

        # If terminal atom belongs to conjugated system, hydrogens are
        # added not to violate rule where atoms should be in one plain.
        my @second_neighbours =
            grep { ! /$atom_id/smx }
            map { @{ $reference_atom_site->{$_}{'connections'} } }
                @connection_ids;

        for my $second_neighbour ( @second_neighbours ) {
            my $second_hybridization =
                $reference_atom_site->{$second_neighbour}{'hybridization'};
            if( $second_hybridization eq 'sp2' ||
                $second_hybridization eq 'sp' ) {
                $side_coord =
                    [ $reference_atom_site->{$second_neighbour}{'Cartn_x'},
                      $reference_atom_site->{$second_neighbour}{'Cartn_y'},
                      $reference_atom_site->{$second_neighbour}{'Cartn_z'}, ];
                last;
            }
        }

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $side_coord,
                                 'global' ) };

        if( scalar @{ $missing_hydrogens } >= 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( -0.5 * $pi ) *
                             sin( 120 * $pi / 180 ) ],
                           [ $bond_length *
                             sin( -0.5 * $pi ) *
                             sin( 120 * $pi / 180 ) ],
                           [ $bond_length *
                             cos( 120 * $pi / 180 ) ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 0.5 * $pi ) *
                             sin( 120 * $pi / 180 ) ],
                           [ $bond_length *
                             sin( 0.5 * $pi ) *
                             sin( 120 * $pi / 180 ) ],
                           [ $bond_length *
                             cos( 120 * $pi / 180 ) ],
                           [ 1 ] ] ] ) };
        }
    }

    return;
}

#
# Adds hydrogens to the atoms which are sp hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp
{
    my ( $parameters, $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;
    my ( $reference_atom_site ) = ( $options->{'reference_atom_site'} );

    $reference_atom_site //= $atom_site;

    my $atom_properties = $parameters->{'_[local]_atom_properties'};
    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $atom_properties->{$atom_type}{'covalent_radius'}{'length'}[1] +
        $atom_properties->{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids = @{ $reference_atom_site->{"$atom_id"}{'connections'} };
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my ( $up_atom_coord,
         $mid_atom_coord,
         $side_coord ) =
             ( $atom_coord{$connection_ids[0]},
               $atom_coord{$atom_id},
               [ $atom_site->{$atom_id}{'Cartn_x'},
                 $atom_site->{$atom_id}{'Cartn_y'} + 1,
                 $atom_site->{$atom_id}{'Cartn_z'}, ], );

    # Generates transformation matrix for transfering atoms to local
    # reference frame.
    my ( $transf_matrix ) =
        @{ switch_ref_frame( $parameters,
                             $mid_atom_coord,
                             $up_atom_coord,
                             $side_coord,
                             'global' ) };

    ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
        @{ matrix_product(
               [ $transf_matrix,
                 [ [ 0 ],
                   [ 0 ],
                   [ - $bond_length ],
                   [ 1 ] ] ] ) };

    return;
}

1;
