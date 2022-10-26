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
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my %options = defined $options ? %{ $options } : ();
    $options{'append_func'} = \&BondProperties::append_rotatable_bonds;
    $options{'mainchain_symbol'} = '';
    $options{'sidechain_symbol'} = 'chi';
    $options{'ignore_terminal_repetition'} = 1;
    my $rotatable_bonds =
        bond_path_search( $parameters, $atom_site, $start_atom_ids, \%options );
    return $rotatable_bonds;
}

sub append_rotatable_bonds
{
    my ( $rotatable_bonds, $atom_site, $atom_id, $parent_atom_ids,
         $include_hetatoms ) = @_;

    my $parent_atom_id = $parent_atom_ids->{$atom_id};
    return if ! defined $parent_atom_id;

    my $is_hetatom = $atom_site->{$atom_id}{'group_PDB'} eq 'HETATM';
    my $is_parent_hetatom =
        $atom_site->{$parent_atom_id}{'group_PDB'} eq 'HETATM';

    # NOTE: this hetatom exception currently will work on single atoms.
    # NOTE: make sure that interaction between 'is_hetatom' and
    # 'include_hetatoms' is correct.
    if( ( ! exists $atom_site->{$atom_id}{'hybridization'} ) &&
        ( ! $is_hetatom ) ) {
        confess "atom with id $atom_id lacks information about " .
            "hybridization"
    }
    if( ( ! exists $atom_site->{$parent_atom_id}{'hybridization'} ) &&
        ( ! $is_parent_hetatom ) ) {
        confess "atom with id $parent_atom_id lacks information about " .
            "hybridization"
    }

    if( $atom_site->{$parent_atom_id}{'hybridization'} eq 'sp3' ||
        $atom_site->{$atom_id}{'hybridization'} eq 'sp3' ||
        ( $include_hetatoms &&
          ( $is_hetatom || $is_parent_hetatom ) &&
          $atom_site->{$atom_id}{'hybridization'} eq '.' ) ) {
        # If last visited atom was sp3, then rotatable bonds from
        # previous atom are copied and the new one is appended.
        push @{ $rotatable_bonds->{$atom_id} },
            [ $parent_atom_id, $atom_id ];
        if( exists $rotatable_bonds->{$parent_atom_id} ) {
            unshift @{ $rotatable_bonds->{$atom_id} },
                @{ $rotatable_bonds->{$parent_atom_id} };
        }
    } else {
        # If last visited atom is sp2 or sp, inherits its rotatable
        # bonds, because double or triple bonds do not rotate.
        if( exists $rotatable_bonds->{$parent_atom_id} ) {
            unshift @{ $rotatable_bonds->{$atom_id} },
                @{ $rotatable_bonds->{$parent_atom_id} };
        }
    }

    return;
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
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my %options = defined $options ? %{ $options } : ();
    $options{'append_func'} = \&BondProperties::append_stretchable_bonds;
    $options{'mainchain_symbol'} = 'd';
    $options{'sidechain_symbol'} = 'r';
    my $stretchable_bonds =
        bond_path_search( $parameters, $atom_site, $start_atom_ids, \%options );
    return $stretchable_bonds;
}

sub append_stretchable_bonds
{
    my ( $bonds, $atom_site, $atom_id, $parent_atom_ids ) = @_;

    my $parent_atom_id = $parent_atom_ids->{$atom_id};
    return if ! defined $parent_atom_id;

    push @{ $bonds->{$atom_id} }, [ $parent_atom_id, $atom_id ];

    # Adds bond if it is a continuation of identified bonds.
    if( exists $bonds->{$parent_atom_id} ) {
        unshift @{ $bonds->{$atom_id} }, @{ $bonds->{$parent_atom_id} };
    }

    return;
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
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my %options = defined $options ? %{ $options } : ();
    $options{'append_func'} = \&BondProperties::append_bendable_angles;
    $options{'mainchain_symbol'} = 'theta';
    $options{'sidechain_symbol'} = 'eta';
    my $bendable_angles =
        bond_path_search( $parameters, $atom_site, $start_atom_ids, \%options );
    return $bendable_angles;
}

sub append_bendable_angles
{
    my ( $angles, $atom_site, $atom_id, $parent_atom_ids ) = @_;

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

sub bond_path_search
{
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my ( $append_func, $ignore_terminal_repetition, $ignore_atoms,
         $include_hetatoms, $ignore_connections ) = (
        $options->{'append_func'},
        $options->{'ignore_terminal_repetition'},
        $options->{'ignore_atoms'},
        $options->{'include_hetatoms'},
        $options->{'ignore_connections'}
    );

    $ignore_terminal_repetition //= 0;
    $ignore_atoms //= {};
    $include_hetatoms //= 0;
    $ignore_connections //= {};

    # By default, N is starting atom for main-chain calculations.
    $start_atom_ids //=
        filter( { 'atom_site' => $atom_site,
                  'include' => { 'label_atom_id' => [ 'N' ] },
                  'data' => [ 'id' ],
                  'is_list' => 1 } );

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

        $append_func->( \%bond_paths, $atom_site, $atom_id, \%parent_atom_ids,
                        $include_hetatoms );

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

            # Depth-first search.
            if( $are_any_sidechain_atoms ) {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
                next;
            }

            # Bread-first search.
            if( $i == 0 ) {
                push @next_atom_ids, $sorted_neighbour_atom_id;
            } else {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
            }
        }
    }

    return name_bond_parameters($parameters, $atom_site, \%bond_paths, $options);
}

sub name_bond_parameters
{
    my ( $parameters, $atom_site, $bonds, $options ) = @_;
    my ( $do_mainchain, $mainchain_symbol, $sidechain_symbol,
         $hetatom_symbol ) = (
        $options->{'do_mainchain'},
        $options->{'mainchain_symbol'},
        $options->{'sidechain_symbol'},
        $options->{'hetatom_symbol'}
    );

    $do_mainchain //= 0;
    $mainchain_symbol //= 'x';
    $sidechain_symbol //= 'y';
    $hetatom_symbol //= '*';

    my %mainchain_atom_names =
        map { $_ => 1 } @{ $parameters->{'_[local]_mainchain_atom_names'} };

    my %bond_parameter_names = ();
    my %visited_bonds = ();
    my %name_counter = (
        "$mainchain_symbol" => 1,
        "$sidechain_symbol" => 1
    );

    for my $atom_ids ( map { @{ $bonds->{$_} } } keys %{ $bonds } ) {
        my $is_visited = $visited_bonds{join(',',@{$atom_ids})};

        next if defined $is_visited && $is_visited;

        $visited_bonds{join(',',@{$atom_ids})} = 1;

        my $are_any_mainchain_atoms =
            any { $mainchain_atom_names{$_} }
            grep { $atom_site->{$_}{'label_atom_id'} }
            @{ $atom_ids };
        my $are_any_hetatoms =
            grep { $atom_site->{$_}{'group_PDB'} eq 'HETATM' } @{ $atom_ids };

        next if ! $do_mainchain && $are_any_mainchain_atoms;

        # Adding parameter names.
        my $bond_parameter_name = "";
        if( $are_any_mainchain_atoms ) {
            $bond_parameter_name .=
                $mainchain_symbol .
                $name_counter{$mainchain_symbol};
            $name_counter{$mainchain_symbol}++;
        }

        if( ! $are_any_mainchain_atoms ) {
            $bond_parameter_name .=
                $sidechain_symbol .
                $name_counter{$sidechain_symbol};
            $name_counter{$sidechain_symbol}++;
        }

        # Adding symbol for hetatoms.
        if( $are_any_hetatoms ) {
            $bond_parameter_name .= $hetatom_symbol;
        }

        $bond_parameter_names{join(',',@{$atom_ids})} =
            $bond_parameter_name;
        $bond_parameter_names{join(',',reverse @{$atom_ids})} =
            $bond_parameter_name;
    }

    my %named_bond_parameters = ();

    for my $atom_id ( keys %{ $bonds } ) {
        for my $atom_ids ( @{ $bonds->{$atom_id} } ) {
            my $bond_parameter_name =
                $bond_parameter_names{join(',',@{$atom_ids})};

            next if ! defined $bond_parameter_name;

            $named_bond_parameters{$atom_id}{$bond_parameter_name} = $atom_ids;
        }
    }

    return \%named_bond_parameters;
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
