package BondProperties;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
    our @EXPORT_OK = qw( bond_type
                         bendable_angles
                         hybridization
                         rotatable_bonds
                         stretchable_bonds
                         unique_rotatables );
}

use Carp;
use List::Util qw( any );
use List::MoreUtils qw( uniq );

use AtomProperties qw( sort_atom_names
                       sort_atom_ids_by_name );
use Measure qw( dihedral_angle );
use PDBxParser qw( filter
                   filter_new );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------- Bond related functions --------------------------- #

#
# Identifies bond type by checking the bond distance.
# Input:
#     ${target,neighbour}_atom - atom data structure (see PDBxParser.pm).
# Output:
#     $bond_type - bond type: single, double, tripple.
#

sub bond_type
{
    my ( $parameters, $target_atom, $neighbour_atom ) = @_;

    my $bond_types = $parameters->{'_[local]_bond_types'};

    my $target_atom_type = $target_atom->{'type_symbol'};
    my $neighbour_atom_type = $neighbour_atom->{'type_symbol'};

    # Precalculates squared distance between atom pairs. Delocalized bonds are
    # described by double or triple bond.
    my $squared_distance =
        ( $neighbour_atom->{'Cartn_x'} - $target_atom->{'Cartn_x'} ) ** 2 +
        ( $neighbour_atom->{'Cartn_y'} - $target_atom->{'Cartn_y'} ) ** 2 +
        ( $neighbour_atom->{'Cartn_z'} - $target_atom->{'Cartn_z'} ) ** 2;

    for my $bond_type ( keys %{ $bond_types } ) {
        if( exists $bond_types->{$bond_type}
                                {$target_atom_type}
                                {$neighbour_atom_type} ) {
            my $bond_length_min = $bond_types->{$bond_type}
                                               {$target_atom_type}
                                               {$neighbour_atom_type}
                                               {'min_length'};
            my $bond_length_max = $bond_types->{$bond_type}
                                               {$target_atom_type}
                                               {$neighbour_atom_type}
                                               {'max_length'};

            if( ( $squared_distance >  $bond_length_min**2 ) &&
                ( $squared_distance <= $bond_length_max**2 ) ) {
                return $bond_type;
            }
        }
    }

    return;
}

#
# Identifies hybridization by examining bond types by distances and the
# hybridizations of the connected atoms.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
#     $options->{'calculate_hybridization'} - ignores hybridizations from
#     parameter file;
# Output:
#     writes 'sp3', 'sp2' or 'sp' value to 'hybridization' key in atom data
#     structure.
#

sub hybridization
{
    # Use connect_atoms before using hybridization function.
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $calculate_hybridization ) = ( $options->{'calculate_hybridization'} );
    $calculate_hybridization //= 0;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $clear_hybridization = $parameters->{'_[local]_clear_hybridization'};

    for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
        # Looks up to a pre-determined hash of known hybridization values and
        # quits early.
        my $atom_name = $atom_site->{$atom_id}{'label_atom_id'};
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        if( exists $clear_hybridization->{$residue_name} &&
            exists $clear_hybridization->{$residue_name}{$atom_name} ) {
            $atom_site->{$atom_id}{'hybridization'} =
                $clear_hybridization->{$residue_name}{$atom_name};
            next;
        }

        # Determines if the connected atoms sits in one plane.
        my $dihedral_angle;
        if( exists $atom_site->{$atom_id}{'connections'} &&
            scalar @{ $atom_site->{$atom_id}{'connections'} } == 3 ) {
            my ( $left_atom_id, $right_atom_id, $up_atom_id ) =
                @{ $atom_site->{$atom_id}{'connections'} };
            $dihedral_angle = dihedral_angle(
                [ map { [ $atom_site->{$_}{'Cartn_x'},
                          $atom_site->{$_}{'Cartn_y'},
                          $atom_site->{$_}{'Cartn_z' } ] }
                  ( $left_atom_id, $right_atom_id, $atom_id,
                    $up_atom_id ) ]
            );

        } elsif( exists $atom_site->{$atom_id}{'connections'} &&
                 scalar @{ $atom_site->{$atom_id}{'connections'} } == 2 ) {
            my ( $left_atom_id, $right_atom_id, $up_atom_id ) =
                @{ $atom_site->{$atom_id}{'connections'} };

            my @neighbours_neighbours;
            for my $neighbours_neighbour (
                ( @{ $atom_site->{$left_atom_id}{'connections'} },
                  @{ $atom_site->{$right_atom_id}{'connections'} } ) ) {
                if( ! any { $neighbours_neighbour eq $_ }
                          ( $left_atom_id, $right_atom_id, $atom_id ) ) {
                    push @neighbours_neighbours, $neighbours_neighbour;
                }
            }
            @neighbours_neighbours = uniq @neighbours_neighbours;

            for my $neighbours_neighbour ( @neighbours_neighbours ) {
                my $current_dihedral_angle = dihedral_angle(
                    [ map { [ $atom_site->{$_}{'Cartn_x'},
                              $atom_site->{$_}{'Cartn_y'},
                              $atom_site->{$_}{'Cartn_z'} ] }
                      ( $left_atom_id, $right_atom_id, $atom_id,
                        $neighbours_neighbour) ]
                );
                if( ( $current_dihedral_angle >=  0.95 * $pi &&
                      $current_dihedral_angle <=  1.05 * $pi ) ||
                    ( $current_dihedral_angle >= -0.05 * $pi &&
                      $current_dihedral_angle <=  0.05 * $pi ) ) {
                    $dihedral_angle = $current_dihedral_angle;
                    last;
                }
            }
        }

        if( defined $dihedral_angle && $dihedral_angle >= 0.95 * $pi ) {
            $atom_site->{$atom_id}{'hybridization'} = 'sp2';
            next;
        }

        # Determines every type of connection.
        my @bond_types;
        for my $connection_id ( @{ $atom_site->{$atom_id}{'connections'} } ) {
            push @bond_types,
                 bond_type( $parameters, $atom_site->{$atom_id},
                            $atom_site->{$connection_id} );
        }

        # Depending on connections, assigns hybridization type.
        if( any { $_ eq 'double' } @bond_types ) {
            $atom_site->{$atom_id}{'hybridization'} = 'sp2';
        } elsif( any { $_ eq 'triple' } @bond_types ) {
            $atom_site->{$atom_id}{'hybridization'} = 'sp';
        } else {
            $atom_site->{$atom_id}{'hybridization'} = 'sp3';
        }
    }

    return;
}

#
# Identifies bonds that can be rotated by torsional angle.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     ${start,next}_atom_id - starting atom id and the next one that is followed
#     in order to identify the direction of the search for rotatable bonds.
# Output:
#     %named_rotatable_bonds - data structure that describes rotatable bonds and
#     the constituent atom ids of the bond. Ex.:
#     {
#       $atom_id_3 => {
#                       $rotatable_bond_name_1 => [ $atom_id_1, $atom_id_2 ],
#                       ...
#                     },
#       ...
#     }
#

sub rotatable_bonds
{
    my ( $atom_site, $start_atom_id, $next_atom_ids, $options ) = @_;
    my ( $ignore_atoms, $include_hetatoms, $ignore_connections ) =
        ( $options->{'ignore_atoms'}, $options->{'include_hetatoms'},
          $options->{'ignore_connections'} );

    $ignore_atoms //= [];
    $include_hetatoms //= 0;
    $ignore_connections //= {};

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' =>[ 'CA' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $next_atom_ids //=  filter( { 'atom_site' => $atom_site,
                                  'include' => { 'label_atom_id' => [ 'CB' ] },
                                  'data' => [ 'id' ],
                                  'is_list' => 1 } );

    if( ! $start_atom_id || ! @{ $next_atom_ids } ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $start_atom_id, @{ $ignore_atoms } );
    my @next_atom_ids = ( @{ $next_atom_ids } );
    my %parent_atom_ids;

    my %rotatable_bonds;

    # Marks parent atom for next atom id.
    for my $next_atom_id ( @{ $next_atom_ids } ) {
        $parent_atom_ids{$next_atom_id} = $start_atom_id;
    }

    # Exists if there are no atoms that is not already visited.
    while( scalar( @next_atom_ids ) != 0 ) {
        # Iterates through every neighbouring atom if it was not visited
        # before.
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my $parent_atom_id = $parent_atom_ids{$atom_id};

            # The direction of bond matters here and is intentional.
            next if $ignore_connections->{$parent_atom_id}{$atom_id};

            my $is_hetatom = $atom_site->{$atom_id}{'group_PDB'} eq 'HETATM';
            my $is_parent_hetatom =
                $atom_site->{$parent_atom_id}{'group_PDB'} eq 'HETATM';

            # NOTE: this hetatom exception currently will work on single atoms.
            # NOTE: make sure that interaction between 'is_hetatom' and
            # 'include_hetatoms' is correct.
            if( ( ! exists $atom_site{$atom_id}{'hybridization'} ) &&
                ( ! $is_hetatom ) ) {
                confess "atom with id $atom_id lacks information about " .
                        "hybridization"
            }
            if( ( ! exists $atom_site{$parent_atom_id}{'hybridization'} ) &&
                ( ! $is_parent_hetatom ) ) {
                confess "atom with id $parent_atom_id lacks information about " .
                        "hybridization"
            }

            if( $atom_site{$parent_atom_id}{'hybridization'} eq 'sp3' ||
                $atom_site{$atom_id}{'hybridization'} eq 'sp3' ||
                ( $include_hetatoms &&
                  ( $is_hetatom || $is_parent_hetatom ) &&
                  $atom_site{$atom_id}{'hybridization'} eq '.' ) ) {
                # If last visited atom was sp3, then rotatable bonds from
                # previous atom are copied and the new one is appended.
                push @{ $rotatable_bonds{$atom_id} },
                     [ $parent_atom_id, $atom_id ];
                if( exists $rotatable_bonds{$parent_atom_id} ) {
                    unshift @{ $rotatable_bonds{$atom_id} },
                            @{ $rotatable_bonds{$parent_atom_id} };
                }
            } else {
                # If last visited atom is sp2 or sp, inherits its rotatable
                # bonds, because double or triple bonds do not rotate.
                if( exists $rotatable_bonds{$parent_atom_id} ) {
                    unshift @{ $rotatable_bonds{$atom_id} },
                            @{ $rotatable_bonds{$parent_atom_id} };
                }
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            if( $include_hetatoms &&
                ( $is_hetatom || $is_parent_hetatom ) &&
                ! exists $atom_site{$atom_id}{'connections_hetatom'} ) {
                confess "atom with id $atom_id lacks 'connections_hetatom' key"
            }

            if( ! $include_hetatoms &&
                ! exists $atom_site{$atom_id}{'connections'} ) {
                confess "atom with id $atom_id lacks 'connections' key"
            }

            # Marks neighbouring atoms.
            if( $include_hetatoms &&
                defined $atom_site{$atom_id}{'connections_hetatom'} ) {
                push @neighbour_atom_ids,
                    @{ $atom_site{$atom_id}{'connections_hetatom'} };
            }
            if( defined $atom_site{$atom_id}{'connections'} ) {
                push @neighbour_atom_ids,
                    @{ $atom_site{$atom_id}{'connections'} };
            }

            # Marks parent atoms for each neighbouring atom.
            for my $neighbour_atom_id ( @neighbour_atom_ids ) {
                if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) &&
                    # HACK: this exception might produce unexpected results.
                    ( ! exists $parent_atom_ids{$neighbour_atom_id} ) &&
                    ( ! $ignore_connections->{$atom_id}{$neighbour_atom_id} ) ) {
                    $parent_atom_ids{$neighbour_atom_id} = $atom_id;
                }
            }
        }

        # Determines next atoms that should be visited.
        @next_atom_ids = (); # Resets value for the new ones to be appended.
        for my $neighbour_atom_id ( uniq @neighbour_atom_ids ) {
            if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) &&
                ( any { $neighbour_atom_id eq $_ } @atom_ids ) ) {
                push @next_atom_ids, $neighbour_atom_id;
            }
        }
    }

    # Removes bonds, if they have the id of the target atom. Also, remove ids,
    # which have no rotatable bonds after previous filtering.
    for my $atom_id ( keys %rotatable_bonds ) {
        my $last_bond_idx = $#{ $rotatable_bonds{$atom_id} };
        if( ( $atom_id == $rotatable_bonds{$atom_id}[$last_bond_idx][0] ||
              $atom_id == $rotatable_bonds{$atom_id}[$last_bond_idx][1] ) ) {
            pop @{ $rotatable_bonds{$atom_id} };
        }
        if( ! @{ $rotatable_bonds{$atom_id} } ) {
            delete $rotatable_bonds{$atom_id};
        }
    }

    # Asigns names for rotatables bonds by first filtering out redundant bonds.
    # TODO: the whole process of naming bonds might be implemented in the while
    # loop above.
    my @unique_bonds;
    for my $bond ( map { @{ $rotatable_bonds{$_} } } keys %rotatable_bonds ) {
        if( ! any { $bond->[0] eq $_->[0] && $bond->[1] eq $_->[1] }
                   @unique_bonds ){
            push @unique_bonds, $bond;
        }
    }

    # Sorts bonds by naming priority.
    my @bond_second_ids = map { $_->[1] } @unique_bonds; # Second atom in the
                                                         # bond.
    my %bond_names;
    my $bond_name_id = 1;
    my $bond_stem = $include_hetatoms ? 'tau' : 'chi';

    if( $include_hetatoms ) {
        # Heteroatom names are sorted according to rotatable bond direction.
        # Names by second atom priority.
        for my $second_atom_id ( @bond_second_ids ) {
            $bond_names{"$second_atom_id"} = "${bond_stem}${bond_name_id}";
            $bond_name_id++;
        }
    } else {
        my @second_names_sorted =
            @{ sort_atom_names(
                   filter( { 'atom_site' => \%atom_site,
                             'include' => { 'id' => \@bond_second_ids },
                             'data' => [ 'label_atom_id' ],
                             'is_list' => 1 } ), { 'sort_type' => 'gn' } ) };

        # Names by second atom priority.
        for my $second_name ( @second_names_sorted ) {
            my $second_atom_id =
                filter( { 'atom_site' => \%atom_site,
                          'include' => { 'label_atom_id' => [ $second_name ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];
            $bond_names{"$second_atom_id"} = "${bond_stem}${bond_name_id}";
            $bond_name_id++;
        }
    }

    # Iterates through rotatable bonds and assigns names by second atom.
    my %named_rotatable_bonds;
    for my $atom_id ( keys %rotatable_bonds ) {
        for my $bond ( @{ $rotatable_bonds{"$atom_id"} } ) {
            my $bond_name = $bond_names{"$bond->[1]"};
                $named_rotatable_bonds{"$atom_id"}{"$bond_name"} = $bond;
        }
    }

    return \%named_rotatable_bonds;
}

#
# Identifies bonds that can be stretched.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     ${start,next}_atom_id - starting atom id and the next one that is followed
#     in order to identify the direction of the search for rotatable bonds.
# Output:
#     %named_stretchable_bonds - data structure that describes stretchable bonds
#     and the constituent atom ids of the bond. Ex.:
#     {
#       $atom_id_3 => {
#                       $stretchable_bond_name_1 => [ $atom_id_1, $atom_id_2 ],
#                       ...
#                     },
#       ...
#     }
#

sub stretchable_bonds
{
    my ( $parameters, $atom_site, $start_atom_id, $options ) = @_;
    my ( $include_hetatoms, $ignore_atoms, $ignore_connections ) =
        ( $options->{'include_hetatoms'}, $options->{'ignore_atoms'},
          $options->{'ignore_connections'} );

    $include_hetatoms //= 0;
    $ignore_atoms //= {};
    $ignore_connections //= {};

    # By default, N is starting atom for main-chain calculations.
    $start_atom_id //=
        filter( { 'atom_site' => $atom_site,
                  'include' => { 'label_atom_id' => [ 'N' ] },
                  'data' => [ 'id' ],
                  'is_list' => 1 } );
    my $stretchable_bonds =
        bond_path_search( $parameters, $atom_site, $start_atom_id,
                          { 'append_func' =>
                                \&BondProperties::append_stretchable_bonds,
                            'naming_func' =>
                                \&BondProperties::name_stretchable_bonds,
                            'include_hetatoms' => $include_hetatoms,
                            'ignore_atoms' => $ignore_atoms,
                            'ignore_connections' => $ignore_connections } );

    return $stretchable_bonds;
}

sub append_stretchable_bonds
{
    my ( $bonds, $atom_id, $parent_atom_ids ) = @_;

    my $parent_atom_id = $parent_atom_ids->{$atom_id};
    return if ! defined $parent_atom_id;

    push @{ $bonds->{$atom_id} }, [ $parent_atom_id, $atom_id ];

    # Adds bond if it is a continuation of identified bonds.
    if( exists $bonds->{$parent_atom_id} ) {
        unshift @{ $bonds->{$atom_id} }, @{ $bonds->{$parent_atom_id} };
    }

    return;
}

sub name_stretchable_bonds
{
    my ( $parameters, $atom_site, $stretchable_bonds, $options ) = @_;
    my ( $do_mainchain, $mainchain_distance_symbol, $sidechain_distance_symbol,
         $hetatom_symbol ) = (
        $options->{'do_mainchain'},
        $options->{'mainchain_distance_symbol'},
        $options->{'sidechain_distance_symbol'},
        $options->{'hetatom_symbol'}
    );

    $do_mainchain //= 0;
    $mainchain_distance_symbol //= 'd';
    $sidechain_distance_symbol //= 'r';
    $hetatom_symbol //= '*';

    my $mainchain_atom_names = $parameters->{'_[local]_mainchain_atom_names'};

    my %bond_names = ();
    my %visited_bonds = ();
    my %bond_counter = (
        "$mainchain_distance_symbol" => 1,
        "$sidechain_distance_symbol" => 1
    );

    for my $atom_ids ( map { @{ $stretchable_bonds->{$_} } }
                       keys %{ $stretchable_bonds } ) {
        my ( $first_atom_id, $second_atom_id ) = @{ $atom_ids };

        next if $visited_bonds{$first_atom_id}{$second_atom_id};

        $visited_bonds{$first_atom_id}{$second_atom_id} = 1;

        my ( $first_atom_name, $second_atom_name ) =
            map { $atom_site->{$_}{'label_atom_id'} }
               @{ $atom_ids };
        my $are_any_mainchain_atoms =
            ( any { $first_atom_name eq $_ } @{ $mainchain_atom_names } ) &&
            ( any { $second_atom_name eq $_ } @{ $mainchain_atom_names } );
        my $are_any_hetatoms =
            grep { $atom_site->{$_}{'group_PDB'} eq 'HETATM' } @{ $atom_ids };

        next if ! $do_mainchain && $are_any_mainchain_atoms;

        # Adding bond names.
        my $bond_name = "";
        if( $are_any_mainchain_atoms ) {
            $bond_name .=
                $mainchain_distance_symbol .
                $bond_counter{$mainchain_distance_symbol};
            $bond_counter{$mainchain_distance_symbol}++;
        }

        if( ! $are_any_mainchain_atoms ) {
            $bond_name .=
                $sidechain_distance_symbol .
                $bond_counter{$sidechain_distance_symbol};
            $bond_counter{$sidechain_distance_symbol}++;
        }

        # Adding symbol for hetatoms.
        if( $are_any_hetatoms ) {
            $bond_name .= $hetatom_symbol;
        }

        $bond_names{$first_atom_id}{$second_atom_id} = $bond_name;
        $bond_names{$second_atom_id}{$first_atom_id} = $bond_name;
    }

    my %named_stretchable_bonds = ();

    for my $atom_id ( keys %{ $stretchable_bonds } ) {
        for my $bond ( @{ $stretchable_bonds->{$atom_id} } ) {
            my $bond_name = $bond_names{$bond->[0]}{$bond->[1]};

            next if ! defined $bond_name;

            $named_stretchable_bonds{$atom_id}{$bond_name} = $bond;
        }
    }

    return \%named_stretchable_bonds;
}

#
# Identifies bonds that can be bent.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     ${start,next}_atom_id - starting atom id and the next one that is followed
#     in order to identify the direction of the search for rotatable bonds.
# Output:
#     %named_bendable_angles - data structure that describes bendable angles
#     and the constituent atom ids of the bond. Ex.:
#     {
#       $atom_id_3 => {
#                       $bendable_angle_name_1 =>
#                           [ $atom_id_1, $atom_id_2, $atom_id_3 ],
#                       ...
#                     },
#       ...
#     }
#

sub bendable_angles
{
    my ( $parameters, $atom_site, $start_atom_id, $options ) = @_;
    my ( $include_hetatoms, $ignore_atoms, $ignore_connections ) =
        ( $options->{'include_hetatoms'}, $options->{'ignore_atoms'},
          $options->{'ignore_connections'} );

    $include_hetatoms //= 0;
    $ignore_atoms //= {};
    $ignore_connections //= {};

    # By default, N is starting atom for main-chain calculations.
    $start_atom_id //=
        filter( { 'atom_site' => $atom_site,
                  'include' => { 'label_atom_id' => [ 'N' ] },
                  'data' => [ 'id' ],
                  'is_list' => 1 } );
    my $bendable_angles =
        bond_path_search( $parameters, $atom_site, $start_atom_id,
                          { 'append_func' =>
                                \&BondProperties::append_bendable_angles,
                            'naming_func' =>
                                \&BondProperties::name_bendable_angles,
                            'include_hetatoms' => $include_hetatoms,
                            'ignore_atoms' => $ignore_atoms,
                            'ignore_connections' => $ignore_connections } );

    return $bendable_angles;
}

sub append_bendable_angles
{
    my ( $angles, $atom_id, $parent_atom_ids ) = @_;

    my $parent_atom_id = $parent_atom_ids->{$atom_id};
    return if ! defined $parent_atom_id;

    my $grandparent_atom_id = $parent_atom_ids->{$parent_atom_id};
    return if ! defined $grandparent_atom_id;

    push @{ $angles->{$atom_id} },
        [ $grandparent_atom_id, $parent_atom_id, $atom_id ];

    # Adds angles if it is a continuation of identified bonds.
    if( exists $angles->{$parent_atom_id} ) {
        unshift @{ $angles->{$atom_id} }, @{ $angles->{$parent_atom_id} };
    }

    return;
}

sub name_bendable_angles
{
    my ( $parameters, $atom_site, $bendable_angles, $options ) = @_;
    my ( $do_mainchain, $mainchain_angle_symbol, $sidechain_angle_symbol,
         $hetatom_symbol ) = (
        $options->{'do_mainchain'},
        $options->{'mainchain_distance_symbol'},
        $options->{'sidechain_distance_symbol'},
        $options->{'hetatom_symbol'}
    );

    $do_mainchain //= 0;
    $mainchain_angle_symbol //= 'theta';
    $sidechain_angle_symbol //= 'eta';
    $hetatom_symbol //= '*';

    my $mainchain_atom_names = $parameters->{'_[local]_mainchain_atom_names'};

    my %angle_names = ();
    my %visited_angles = ();
    my %angle_counter = (
        "$mainchain_angle_symbol" => 1,
        "$sidechain_angle_symbol" => 1
    );

    for my $atom_ids ( map { @{ $bendable_angles->{$_} } }
                       keys %{ $bendable_angles } ) {
        my ( $first_atom_id, $second_atom_id, $third_atom_id ) = @{ $atom_ids };

        next if $visited_angles{$first_atom_id}{$second_atom_id}{$third_atom_id};
        $visited_angles{$first_atom_id}{$second_atom_id}{$third_atom_id} = 1;

        my ( $first_atom_name, $second_atom_name, $third_atom_name ) =
            map { $atom_site->{$_}{'label_atom_id'} }
               @{ $atom_ids };
        my $are_any_mainchain_atoms =
            ( any { $first_atom_name eq $_ }  @{ $mainchain_atom_names } ) &&
            ( any { $second_atom_name eq $_ } @{ $mainchain_atom_names } ) &&
            ( any { $third_atom_name eq $_ }  @{ $mainchain_atom_names } );
        my $are_any_hetatoms =
            grep { $atom_site->{$_}{'group_PDB'} eq 'HETATM' } @{ $atom_ids };

        next if ! $do_mainchain && $are_any_mainchain_atoms;

        # Adding angle names.
        my $angle_name = "";
        if( $are_any_mainchain_atoms ) {
            $angle_name .=
                $mainchain_angle_symbol .
                $angle_counter{$mainchain_angle_symbol};
            $angle_counter{$mainchain_angle_symbol}++;
        }

        if( ! $are_any_mainchain_atoms ) {
            $angle_name .=
                $sidechain_angle_symbol .
                $angle_counter{$sidechain_angle_symbol};
            $angle_counter{$sidechain_angle_symbol}++;
        }

        # Adding symbol for hetatoms.
        if( $are_any_hetatoms ) {
            $angle_name .= $hetatom_symbol;
        }

        $angle_names{$first_atom_id}{$second_atom_id}{$third_atom_id} =
            $angle_name;
        $angle_names{$third_atom_id}{$second_atom_id}{$first_atom_id} =
            $angle_name;
    }

    my %named_bendable_angles = ();

    for my $atom_id ( keys %{ $bendable_angles } ) {
        for my $angle ( @{ $bendable_angles->{$atom_id} } ) {
            my $angle_name = $angle_names{$angle->[0]}{$angle->[1]}{$angle->[2]};

            next if ! defined $angle_name;

            $named_bendable_angles{$atom_id}{$angle_name} = $angle;
        }
    }

    return \%named_bendable_angles;
}

sub bond_path_search
{
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my ( $append_func, $naming_func, $ignore_atoms, $include_hetatoms,
         $ignore_connections ) = (
        $options->{'append_func'},
        $options->{'naming_func'},
        $options->{'ignore_atoms'},
        $options->{'include_hetatoms'},
        $options->{'ignore_connections'}
    );

    $ignore_atoms //= {};
    $include_hetatoms //= 0;
    $ignore_connections //= {};

    if( ! @{ $start_atom_ids } ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my %visited_atom_ids = %{ $ignore_atoms };
    my @next_atom_ids = ( @{ $start_atom_ids } );
    my %parent_atom_ids;

    my $mainchain_atom_names = $parameters->{'_[local]_mainchain_atom_names'};

    my %bond_paths = ();

    # Exists if there are no atoms that is not already visited.
    while( @next_atom_ids ) {
        my ( $atom_id ) = pop @next_atom_ids;

        next if $visited_atom_ids{$atom_id};
        $visited_atom_ids{$atom_id} = 1;

        $append_func->( \%bond_paths, $atom_id, \%parent_atom_ids );

        # Marks neighbouring atoms.
        my @neighbour_atom_ids = ();
        if( defined $atom_site{$atom_id}{'connections'} ) {
            push @neighbour_atom_ids,
                grep { ! $ignore_connections->{$atom_id}{$_} }
                grep { ! $ignore_atoms->{$_} }
                @{ $atom_site{$atom_id}{'connections'} };
        }
        if( $include_hetatoms &&
            defined $atom_site{$atom_id}{'connections_hetatom'} ) {
            push @neighbour_atom_ids,
                grep { ! $ignore_connections->{$atom_id}{$_} }
                grep { ! $ignore_atoms->{$_} }
                @{ $atom_site{$atom_id}{'connections_hetatom'} };
        }

        my @sorted_neighbour_atom_ids =
            @{ sort_atom_ids_by_name( \@neighbour_atom_ids, \%atom_site ) };

        for( my $i = 0; $i <= $#sorted_neighbour_atom_ids; $i++ ) {
            my $sorted_neighbour_atom_id = $sorted_neighbour_atom_ids[$i];

            next if $ignore_atoms->{$sorted_neighbour_atom_id};
            next if $ignore_connections->{$atom_id}{$sorted_neighbour_atom_id};
            next if $visited_atom_ids{$sorted_neighbour_atom_id};

            $parent_atom_ids{$sorted_neighbour_atom_id} = $atom_id;

            # Depending on if it is mainchain or sidechain bonds, the bond
            # search changes from deapth-first search to breadth-first search
            # accordingly.
            my $are_any_sidechain_atoms =
                ! ( any { $atom_site->{$sorted_neighbour_atom_id}
                                      {'label_atom_id'} eq $_ }
                       @{ $mainchain_atom_names } ) ||
                ! ( any { $atom_site->{$sorted_neighbour_atom_id}
                                      {'label_atom_id'} eq $_ }
                       @{ $mainchain_atom_names } );

            # Bread-first search.
            if( $are_any_sidechain_atoms ) {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
                next;
            }

            if( $i == 0 ) {
                push @next_atom_ids, $sorted_neighbour_atom_id;
            } else {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
            }
        }
    }

    return $naming_func->( $parameters, $atom_site, \%bond_paths );
}

#
# Identifies unique rotatable bonds in selected group of residue atoms.
# Input:
#     $atom_site - atom data structure (see PDBxParser.pm).
# Output:
#     %unique_rotatable_bonds - hash of arrays that point to unique rotatable
#     bonds.
#     Ex.: { 'chi1' => [ 1, 2 ],
#            'chi2' => [ 2, 3 ] }
#

sub unique_rotatables
{
    my ( $atom_site ) = @_;

    my $rotatable_bonds = rotatable_bonds( $atom_site );

    my %unique_rotatable_bonds;
    for my $atom_id ( keys %{ $rotatable_bonds } ) {
        for my $angle_name ( keys %{ $rotatable_bonds->{"$atom_id"} } ){
            if( ! exists $unique_rotatable_bonds{"$angle_name"} ) {
                $unique_rotatable_bonds{"$angle_name"} =
                    $rotatable_bonds->{"$atom_id"}{"$angle_name"};
            }
        }
    }

    return \%unique_rotatable_bonds;
}

1;
