package BondProperties;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
    our @EXPORT_OK = qw( bond_type
                         %COVALENT_BOND_COMB
                         covalent_bond_combinations
                         hybridization
                         rotatable_bonds
                         unique_rotatables );
}

use List::Util qw( any
                   min );
use List::MoreUtils qw( uniq );

use AtomProperties qw( %ATOMS
                       sort_atom_names );
use ConnectAtoms qw( connect_atoms );
use Constants qw( $PI );
use Measure qw( dihedral_angle );
use MoleculeProperties qw( %CLEAR_HYBRIDIZATION );
use PDBxParser qw( filter );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------------- Constants ---------------------------------- #

my %BOND_TYPES = (
    'single' => {
        'H' => {
            'H' => {
                'min_length' =>
                    2 * ( $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] -
                          $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] ),
                'max_length' =>
                    2 * ( $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                          $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] ),
            },
            'C' => {
                'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0],
            },
            'N' => {
                'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
            },
            'O' => {
                'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
            },
            'S' => {
                'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
            },
        },
        'C' => {
            'C' => {
                'min_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] -
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] ),
                'max_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] ),
            },
            'N' => {
                'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
            },
            'O' => {
                'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
            },
            'S' => {
                'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
            },
        },
        'N' => {
            'N' => {
                'min_length' =>
                    2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] -
                          $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] ),
                'max_length' =>
                    2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                          $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] ),
            },
            'O' => {
                'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
            },
            'S' => {
                'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
            },
        },
        'O' => {
            'O' => {
                'min_length' =>
                    2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] -
                          $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] ),
                'max_length' =>
                    2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                          $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] ),
            },
            'S' => {
                'min_length' => $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] -
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
                'max_length' => $ATOMS{'O'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] +
                                $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
            },
        },
        'S' => {
            'S' => {
                'min_length' =>
                    2 * ( $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] -
                          $ATOMS{'S'}{'covalent_radius'}{'error'}->[0] ),
                'max_length' =>
                    2 * ( $ATOMS{'S'}{'covalent_radius'}{'length'}->[0] +
                          $ATOMS{'S'}{'covalent_radius'}{'error'}->[0] ),
            },
        },
    },
    'double' => {
        'C' => {
            'C' => {
                'min_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] -
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] ),
                'max_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] +
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] ),
            },
            'N' => {
                'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[1],
                'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[1],
            },
            'O' => {
                'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] -
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
                'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
            },
        },
        'N' => {
            'N' => {
                'min_length' =>
                    2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] -
                          $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] ),
                'max_length' =>
                    2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] +
                          $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] ),
            },
            'O' => {
                'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] -
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] -
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
                'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] +
                                $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] +
                                $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
            },
        },
        'O' => {
            'O' => {
                'min_length' =>
                    2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] -
                          $ATOMS{'O'}{'covalent_radius'}{'error'}->[1] ),
                'max_length' =>
                    2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[1] +
                          $ATOMS{'O'}{'covalent_radius'}{'error'}->[1] ),
            },
        },
    },
    'triple' => {
        'C' => {
            'C' => {
                'min_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[2] -
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[2] ),
                'max_length' =>
                    2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[2] +
                          $ATOMS{'C'}{'covalent_radius'}{'error'}->[2] ),
            },
        },
    },
);

our %COVALENT_BOND_COMB = %{ covalent_bond_combinations( \%ATOMS ) };

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
    my ( $target_atom, $neighbour_atom ) = @_;

    my $target_atom_type = $target_atom->{'type_symbol'};
    my $neighbour_atom_type = $neighbour_atom->{'type_symbol'};

    # Precalculates squared distance between atom pairs. Delocalized bonds are
    # described by double or triple bond.
    my $squared_distance =
        ( $neighbour_atom->{'Cartn_x'} - $target_atom->{'Cartn_x'} ) ** 2 +
        ( $neighbour_atom->{'Cartn_y'} - $target_atom->{'Cartn_y'} ) ** 2 +
        ( $neighbour_atom->{'Cartn_z'} - $target_atom->{'Cartn_z'} ) ** 2;

    for my $bond_type ( keys %BOND_TYPES ) {
        if( exists $BOND_TYPES{$bond_type}
                              {$target_atom_type}
                              {$neighbour_atom_type} ||
            exists $BOND_TYPES{$bond_type}
                              {$neighbour_atom_type}
                              {$target_atom_type} ) {
            my $bond_length_min = $BOND_TYPES{$bond_type}
                                             {$target_atom_type}
                                             {$neighbour_atom_type}
                                             {'min_length'} ||
                                  $BOND_TYPES{$bond_type}
                                             {$neighbour_atom_type}
                                             {$target_atom_type}
                                             {'min_length'};
            my $bond_length_max = $BOND_TYPES{$bond_type}
                                             {$target_atom_type}
                                             {$neighbour_atom_type}
                                             {'max_length'} ||
                                  $BOND_TYPES{$bond_type}
                                             {$neighbour_atom_type}
                                             {$target_atom_type}
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
# Determines all covalent bond combinations.
# Inputs:
#     \%ATOMS - atom parameter data structure.
# Outputs:
#     %covalent_bond_combinations - data structure for storing covalent bond
#     combinations.
#     Ex.: { 'C' =>
#                { 'C' => { 'length' => [ 0.1, 0.2 ],
#                           'error'  => [ 0.01, 0.01 ] } } }
#

sub covalent_bond_combinations
{
    my ( $ATOMS ) = @_;

    my %covalent_bond_combinations;
    for my $atom_i ( keys %ATOMS ) {
        for my $atom_j ( keys %ATOMS ) {
            for my $i ( 0..min( $#{ $ATOMS{$atom_i}{'covalent_radius'}{'length'} },
                                $#{ $ATOMS{$atom_j}{'covalent_radius'}{'length'} } ) ) {
                push @{ $covalent_bond_combinations{$atom_i}{$atom_j}{'length'} },
                    $ATOMS{$atom_i}{'covalent_radius'}{'length'}[$i] +
                    $ATOMS{$atom_j}{'covalent_radius'}{'length'}[$i];
                push @{ $covalent_bond_combinations{$atom_i}{$atom_j}{'error'} },
                    $ATOMS{$atom_i}{'covalent_radius'}{'error'}[$i] +
                    $ATOMS{$atom_j}{'covalent_radius'}{'error'}[$i];
            }
        }
    }

    return \%covalent_bond_combinations;
}

#
# Identifies hybridization by examining bond types by distances and the
# hybridizations of the connected atoms.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
# Output:
#     writes 'sp3', 'sp2' or 'sp' value to 'hybridization' key in atom data
#     structure.
#

sub hybridization
{
    # Use connect_atoms before using hybridization function.
    my ( $atom_site ) = @_;

    for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
        # Looks up to a pre-determined hash of known hybridization values and
        # quits early.
        my $atom_name = $atom_site->{$atom_id}{'label_atom_id'};
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        if( exists $CLEAR_HYBRIDIZATION{$residue_name} &&
            exists $CLEAR_HYBRIDIZATION{$residue_name}{$atom_name} ) {
            $atom_site->{$atom_id}{'hybridization'} =
                $CLEAR_HYBRIDIZATION{$residue_name}{$atom_name};
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
                if( ( $current_dihedral_angle >=  0.95 * $PI &&
                      $current_dihedral_angle <=  1.05 * $PI ) ||
                    ( $current_dihedral_angle >= -0.05 * $PI &&
                      $current_dihedral_angle <=  0.05 * $PI ) ) {
                    $dihedral_angle = $current_dihedral_angle;
                    last;
                }
            }
        }

        if( defined $dihedral_angle && $dihedral_angle >= 0.95 * $PI ) {
            $atom_site->{$atom_id}{'hybridization'} = 'sp2';
            next;
        }

        # Determines every type of connection.
        my @bond_types;
        for my $connection_id ( @{ $atom_site->{$atom_id}{'connections'} } ) {
            push @bond_types,
                 bond_type( $atom_site->{$atom_id},
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
    my ( $atom_site, $start_atom_id, $next_atom_id ) = @_;

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' => [ 'CA' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];
    $next_atom_id //=  filter( { 'atom_site' => $atom_site,
                                 'include' => { 'label_atom_id' => [ 'CB' ] },
                                 'data' => [ 'id' ],
                                 'is_list' => 1 } )->[0];

    if( ! $start_atom_id || ! $next_atom_id ) { return {}; }

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $start_atom_id );
    my @next_atom_ids = ( $next_atom_id );
    my %parent_atom_ids;

    my %rotatable_bonds;

    # Marks parent atom for next atom id.
    $parent_atom_ids{$next_atom_id} = $start_atom_id;

    # Connects and determines hybridization for each atom.
    connect_atoms( \%atom_site );
    hybridization( \%atom_site );

    # Exists if there are no atoms that is not already visited.
    while( scalar( @next_atom_ids ) != 0 ) {
        # Iterates through every neighbouring atom if it was not visited
        # before.
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my $parent_atom_id = $parent_atom_ids{$atom_id};

            if( $atom_site{$parent_atom_id}{'hybridization'} eq 'sp3' ||
                $atom_site{$atom_id}{'hybridization'} eq 'sp3' ) {
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

            # Marks neighbouring atoms.
            push @neighbour_atom_ids,
                 @{ $atom_site{$atom_id}{'connections'} };

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
