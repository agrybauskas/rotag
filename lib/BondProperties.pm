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

use AtomProperties qw( sort_atom_names );
use Measure qw( dihedral_angle );
use PDBxParser qw( filter );
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
            $dihedral_angle =
                    dihedral_angle(
                        [ [ $atom_site->{$left_atom_id}{'Cartn_x'},
                            $atom_site->{$left_atom_id}{'Cartn_y'},
                            $atom_site->{$left_atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$right_atom_id}{'Cartn_x'},
                            $atom_site->{$right_atom_id}{'Cartn_y'},
                            $atom_site->{$right_atom_id}{'Cartn_z'}],
                          [ $atom_site->{$atom_id}{'Cartn_x'},
                            $atom_site->{$atom_id}{'Cartn_y'},
                            $atom_site->{$atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$up_atom_id}{'Cartn_x'},
                            $atom_site->{$up_atom_id}{'Cartn_y'},
                            $atom_site->{$up_atom_id}{'Cartn_z'} ],
                        ] );

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
                my $current_dihedral_angle =
                    dihedral_angle(
                        [ [ $atom_site->{$left_atom_id}{'Cartn_x'},
                            $atom_site->{$left_atom_id}{'Cartn_y'},
                            $atom_site->{$left_atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$right_atom_id}{'Cartn_x'},
                            $atom_site->{$right_atom_id}{'Cartn_y'},
                            $atom_site->{$right_atom_id}{'Cartn_z'}],
                          [ $atom_site->{$atom_id}{'Cartn_x'},
                            $atom_site->{$atom_id}{'Cartn_y'},
                            $atom_site->{$atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$neighbours_neighbour}{'Cartn_x'},
                            $atom_site->{$neighbours_neighbour}{'Cartn_y'},
                            $atom_site->{$neighbours_neighbour}{'Cartn_z'} ],
                        ] );
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
    my ( $atom_site, $start_atom_id, $next_atom_id, $options ) = @_;
    my ( $do_hetatoms, $ignore_connections ) =
        ( $options->{'do_hetatoms'}, $options->{'ignore_connections'} );

    $do_hetatoms //= 0;
    $ignore_connections = [];

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { 'atom_site' => $atom_site,
                                 'include' =>
                                 { $do_hetatoms ?
                                   ( 'label_atom_id' =>[ 'N' ] ) :
                                   ( 'label_atom_id' =>[ 'CA' ] ) },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $next_atom_id //=  filter( { 'atom_site' => $atom_site,
                                 'include' =>
                                 { $do_hetatoms ?
                                   ( 'label_atom_id' => [ 'CA' ] ) :
                                   ( 'label_atom_id' => [ 'CB' ] ) },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];

    if( ! $start_atom_id || ! $next_atom_id ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $start_atom_id, @{ $ignore_connections } );
    my @next_atom_ids = ( $next_atom_id );
    my %parent_atom_ids;

    my %rotatable_bonds;

    # Marks parent atom for next atom id.
    $parent_atom_ids{$next_atom_id} = $start_atom_id;

    # Exists if there are no atoms that is not already visited.
    while( scalar( @next_atom_ids ) != 0 ) {
        # Iterates through every neighbouring atom if it was not visited
        # before.
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my $parent_atom_id = $parent_atom_ids{$atom_id};

            if( ( ! exists $atom_site{$atom_id}{'hybridization'} ) &&
                ( ! $do_hetatoms ) ) {
                confess "atom with id $atom_id lacks information about " .
                        "hybridization"
            }
            if( ( ! exists $atom_site{$parent_atom_id}{'hybridization'} ) &&
                ( ! $do_hetatoms ) ) {
                confess "atom with id $parent_atom_id lacks information about " .
                        "hybridization"
            }

            if( $atom_site{$parent_atom_id}{'hybridization'} eq 'sp3' ||
                $atom_site{$atom_id}{'hybridization'} eq 'sp3' ||
                ($atom_site{$atom_id}{'hybridization'} eq '.' && $do_hetatoms) ){
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

            if( ! exists $atom_site{$atom_id}{'connections'} ) {
                confess "atom with id $atom_id lacks 'connections' key"
            }

            # Marks neighbouring atoms.
            push @neighbour_atom_ids, @{ $atom_site{$atom_id}{'connections'} };

            # Marks parent atoms for each neighbouring atom.
            for my $neighbour_atom_id ( @neighbour_atom_ids ) {
                if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) &&
                    # HACK: this exception might produce unexpected results.
                    ( ! exists $parent_atom_ids{$neighbour_atom_id} ) ) {
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
    my @second_names_sorted =
        @{ sort_atom_names(
               filter( { 'atom_site' => \%atom_site,
                         'include' => { 'id' => \@bond_second_ids },
                         'data' => [ 'label_atom_id' ],
                         'is_list' => 1 } ), { 'sort_type' => 'gn' } ) };

    my %bond_names; # Names by second atom priority.
    my $bond_name_id = 1;
    for my $second_name ( @second_names_sorted ) {
        my $second_atom_id =
            filter( { 'atom_site' => \%atom_site,
                      'include' => { 'label_atom_id' => [ $second_name ] },
                      'data' => [ 'id' ],
                      'is_list' => 1 } )->[0];
        $bond_names{"$second_atom_id"} = "chi$bond_name_id";
        $bond_name_id++;
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
    # TODO: it seems that $next_atom_ids should be applied to bond rotation and
    # angle bending functions.
    my ( $atom_site, $start_atom_id, $next_atom_ids ) = @_;

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' => [ 'CA' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $next_atom_ids //=  filter( { 'atom_site' => $atom_site,
                                  'include' => { 'label_atom_id' => [ 'CB' ] },
                                  'data' => [ 'id' ],
                                  'is_list' => 1 } );

    if( ! $start_atom_id || ! @{ $next_atom_ids } ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $start_atom_id );
    my @next_atom_ids = ( @{ $next_atom_ids } );
    my %parent_atom_ids;

    my %stretchable_bonds;

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

            push @{ $stretchable_bonds{$atom_id} },
                [ $parent_atom_id, $atom_id ];
            if( exists $stretchable_bonds{$parent_atom_id} ) {
                unshift @{ $stretchable_bonds{$atom_id} },
                    @{ $stretchable_bonds{$parent_atom_id} };
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            if( ! exists $atom_site{$atom_id}{'connections'} ) {
                confess "atom with id $atom_id lacks 'connections' key"
            }

            # Marks neighbouring atoms.
            push @neighbour_atom_ids, @{ $atom_site{$atom_id}{'connections'} };

            # Marks parent atoms for each neighbouring atom.
            for my $neighbour_atom_id ( @neighbour_atom_ids ) {
                if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) &&
                    # HACK: this exception might produce unexpected results.
                    ( ! exists $parent_atom_ids{$neighbour_atom_id} ) ) {
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

    # Asigns names for stretchable bonds by first filtering out redundant bonds.
    # TODO: the whole process of naming bonds might be implemented in the while
    # loop above.
    my @unique_bonds;
    for my $bond ( map { @{ $stretchable_bonds{$_} } } keys %stretchable_bonds ) {
        if( ! any { $bond->[0] eq $_->[0] && $bond->[1] eq $_->[1] }
                   @unique_bonds ){
            push @unique_bonds, $bond;
        }
    }

    # Sorts bonds by naming priority.
    my @bond_second_ids = map { $_->[1] } @unique_bonds; # Second atom in the
                                                         # bond.
    my @second_names_sorted =
        @{ sort_atom_names(
               filter( { 'atom_site' => \%atom_site,
                         'include' => { 'id' => \@bond_second_ids },
                         'data' => [ 'label_atom_id' ],
                         'is_list' => 1 } ), { 'sort_type' => 'gn' } ) };

    my %bond_names; # Names by second atom priority.
    my $bond_name_id = 1;
    for my $second_name ( @second_names_sorted ) {
        my $second_atom_id =
            filter( { 'atom_site' => \%atom_site,
                      'include' => { 'label_atom_id' => [ $second_name ] },
                      'data' => [ 'id' ],
                      'is_list' => 1 } )->[0];
        $bond_names{"$second_atom_id"} = "r$bond_name_id";
        $bond_name_id++;
    }

    # Iterates through stretchable bonds and assigns names by second atom.
    my %named_stretchable_bonds;
    for my $atom_id ( keys %stretchable_bonds ) {
        for my $bond ( @{ $stretchable_bonds{"$atom_id"} } ) {
            my $bond_name = $bond_names{"$bond->[1]"};
            $named_stretchable_bonds{"$atom_id"}{"$bond_name"} = $bond;
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
    my ( $atom_site, $start_atom_id, $next_atom_id, $previous_atom_id ) = @_;

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' => [ 'CA' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $next_atom_id //=  filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' => [ 'CB' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $previous_atom_id //=  filter( { 'atom_site' => $atom_site,
                                     'include' => { 'label_atom_id' => [ 'N' ] },
                                     'data' => [ 'id' ],
                                     'is_list' => 1 } )->[0];

    if( ! $start_atom_id || ! $next_atom_id ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $previous_atom_id, $start_atom_id );
    my @next_atom_ids = ( $next_atom_id );
    my %parent_atom_ids;

    my %bendable_angles;

    # Marks parent and grandparent atoms for next atom id.
    $parent_atom_ids{$next_atom_id} = $start_atom_id;
    $parent_atom_ids{$start_atom_id} = $previous_atom_id;

    # Exists if there are no atoms that is not already visited.
    while( scalar( @next_atom_ids ) != 0 ) {
        # Iterates through every neighbouring atom if it was not visited
        # before.
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my $parent_atom_id = $parent_atom_ids{$atom_id};

            # Detects grandparent atom if it exists.
            my $grandparent_atom_id;
            if( exists $parent_atom_ids{$parent_atom_id} ) {
                $grandparent_atom_id = $parent_atom_ids{$parent_atom_id};
            }

            push @{ $bendable_angles{$atom_id} },
                [ $grandparent_atom_id, $parent_atom_id, $atom_id ];
            if( exists $bendable_angles{$parent_atom_id} ) {
                unshift @{ $bendable_angles{$atom_id} },
                    @{ $bendable_angles{$parent_atom_id} };
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            if( ! exists $atom_site{$atom_id}{'connections'} ) {
                confess "atom with id $atom_id lacks 'connections' key"
            }

            # Marks neighbouring atoms.
            push @neighbour_atom_ids, @{ $atom_site{$atom_id}{'connections'} };

            # Marks parent atoms for each neighbouring atom.
            for my $neighbour_atom_id ( @neighbour_atom_ids ) {
                if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) &&
                    # HACK: this exception might produce unexpected results.
                    ( ! exists $parent_atom_ids{$neighbour_atom_id} ) ) {
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

    # Asigns names for bendable angles by first filtering out redundant angles.
    # TODO: the whole process of naming bonds might be implemented in the while
    # loop above.
    my @unique_bonds;
    for my $bond ( map { @{ $bendable_angles{$_} } } keys %bendable_angles ) {
        if( ! any { $bond->[0] eq $_->[0] && $bond->[1] eq $_->[1] }
                   @unique_bonds ){
            push @unique_bonds, $bond;
        }
    }

    # Sorts bonds by naming priority.
    my @bond_second_ids = map { $_->[1] } @unique_bonds; # Second atom in the
                                                         # bond.
    my @bond_third_ids  = map { $_->[2] } @unique_bonds; # Second atom in the
                                                         # bond.

    my @second_names_sorted =
        @{ sort_atom_names(
               filter( { 'atom_site' => \%atom_site,
                         'include' => { 'id' => \@bond_second_ids },
                         'data' => [ 'label_atom_id' ],
                         'is_list' => 1 } ), { 'sort_type' => 'gn' } ) };
    my @third_names_sorted =
        @{ sort_atom_names(
               filter( { 'atom_site' => \%atom_site,
                         'include' => { 'id' => \@bond_third_ids },
                         'data' => [ 'label_atom_id' ],
                         'is_list' => 1 } ), { 'sort_type' => 'gn' } ) };

    my %angle_names; # Names by second atom priority.
    my $angle_name_id = 1;
    for my $second_name ( @second_names_sorted ) {
        my $second_atom_id =
            filter( { 'atom_site' => \%atom_site,
                      'include' => { 'label_atom_id' => [ $second_name ] },
                      'data' => [ 'id' ],
                      'is_list' => 1 } )->[0];
        $angle_names{"$second_atom_id"} = "theta$angle_name_id";
        $angle_name_id++;
    }

    # Iterates through bendable angles and assigns names by second atom.
    my %named_bendable_angles;
    for my $atom_id ( keys %bendable_angles ) {
        for my $angle ( @{ $bendable_angles{"$atom_id"} } ) {
            my $angle_name = $angle_names{"$angle->[1]"};
            $named_bendable_angles{"$atom_id"}{"$angle_name"} = $angle;
        }
    }

    return \%named_bendable_angles;
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
