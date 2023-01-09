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
use BondPath;
use Measure qw( dihedral_angle );
use PDBxParser qw( filter_new );
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
#     $start_atom_ids - starting atom ids;
#     $options->{'include_hetatoms'} - includes heteroatoms.
# Output:
#     %rotatable_bonds - data structure that describes rotatable bonds and
#     the constituent atom ids of the bond. Ex.:
#     {
#       $atom_id_3 => {
#                       $rotatable_bond_name_1 => [
#                           $atom_id_1, $atom_id_2, $atom_id_3, $atom_id_4
#                       ],
#                       ...
#                     },
#       ...
#     }
#

sub rotatable_bonds
{
    my ( $parameters, $atom_site, $start_atom_ids, $options ) = @_;
    my ( $bond_paths, $include_mainchain, $include_hetatoms ) = (
        $options->{'bond_paths'},
        $options->{'include_mainchain'},
        $options->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $explicit_dihedral_names = $parameters->{'_[local]_dihedral_angle_name'};
    # TODO: maybe it should be moved to force field parameter file?
    my $ignore_connections = {
        'label_atom_id' => {
            'N' => { 'CD' => 1, # For PRO.
                     'C' => 1 },
            'C' => { 'O' => 1,
                     ( $include_hetatoms ? ( 'CA' => 1 ) : ( 'N' => 1 ) ) },
        },
    };

    if( $include_hetatoms ) {
        $start_atom_ids = filter_new(
            $atom_site,
            { 'include' => { 'label_atom_id' => [ 'C' ] },
              'return_data' => 'id'
        } );
    }

    $bond_paths //= BondPath->new( {
        'atom_site' => $atom_site,
        'start_atom_ids' => $start_atom_ids,
        'include_hetatoms' => $include_hetatoms,
        'ignore_connections' => $ignore_connections,
    } );

    my %rotatable_bonds = ();
    my %parent_atom_ids = ();
    my %shared_bonds = ();
    my %order = ();
    for my $order ( sort { $a <=> $b } keys %{ $bond_paths } ) {
        my ( $third_atom_id ) = keys %{ $bond_paths->{$order} };
        my $fourth_atom_id = $bond_paths->{$order}{$third_atom_id};

        $parent_atom_ids{$fourth_atom_id} = $third_atom_id;
        $order{$third_atom_id}{$fourth_atom_id} = $order;

        my $second_atom_id = $parent_atom_ids{$third_atom_id};

        next if ! defined $second_atom_id;

        my $first_atom_id = $parent_atom_ids{$second_atom_id};

        next if ! defined $first_atom_id;

        # Checks for mainchains and heteroatoms.
        next if ! $include_mainchain &&
            ! contains_sidechain_atoms( $parameters,
                                        $atom_site,
                                        [ $first_atom_id, $second_atom_id,
                                          $third_atom_id, $fourth_atom_id ] ) &&
            ! contains_hetatoms( $atom_site,
                                 [ $first_atom_id, $second_atom_id,
                                   $third_atom_id, $fourth_atom_id ] );

        # Check on hybridization.
        if( ! exists $atom_site->{$second_atom_id}{'hybridization'} ) {
            confess "atom with id $second_atom_id lacks information about " .
                "hybridization";
        }
        if( ! exists $atom_site->{$third_atom_id}{'hybridization'} ){
            confess "atom with id $third_atom_id lacks information " .
                "about hybridization";
        }

        # Appending dihedral angles.
        if( $atom_site->{$second_atom_id}{'hybridization'} eq 'sp3' ||
            $atom_site->{$third_atom_id}{'hybridization'} eq 'sp3' ||
            ( $include_hetatoms &&
              $atom_site->{$fourth_atom_id}{'group_PDB'} eq 'HETATM' &&
              $atom_site->{$fourth_atom_id}{'hybridization'} eq '.' ) ) {
            # If at one of the bond atoms are sp3, it is rotatable.
            if( exists $shared_bonds{$second_atom_id}{$third_atom_id} ) {
                push @{ $rotatable_bonds{$fourth_atom_id} },
                    $shared_bonds{$second_atom_id}{$third_atom_id};
            } else {
                push @{ $rotatable_bonds{$fourth_atom_id} },
                    [ $first_atom_id, $second_atom_id, $third_atom_id,
                      $fourth_atom_id ];
            }
        }

        # If bond atoms are sp2/sp (do not rotate) or just is a continuation of
        # the bond chain, inherits its previous atom's rotatable bonds.
        if( exists $rotatable_bonds{$third_atom_id} ) {
            unshift @{ $rotatable_bonds{$fourth_atom_id} },
                @{ $rotatable_bonds{$third_atom_id} };
        }

        # Specific dihedral angle are unique and are described by one set of
        # dihedral angles that are determined by the order (name hierarchy).
        if( ! exists $shared_bonds{$second_atom_id}{$third_atom_id} ) {
            $shared_bonds{$second_atom_id}{$third_atom_id} = [
                $first_atom_id, $second_atom_id, $third_atom_id, $fourth_atom_id
            ];
        }
    }

    # Naming the rotatable bonds.
    my %named_rotatable_bonds = ();
    for my $atom_id ( keys %rotatable_bonds ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $rotatable_bonds{$atom_id} } ) {
            my $rotatable_bond_name =
                join '-', map { $atom_site->{$_}{'label_atom_id'} }
                             @{ $bond_atom_ids };
            my $simplified_rotatable_bond_name =
                join '-', ( '.',
                            $atom_site->{$bond_atom_ids->[1]}{'label_atom_id'},
                            $atom_site->{$bond_atom_ids->[2]}{'label_atom_id'},
                            '.');

            if( exists $explicit_dihedral_names->{$residue_name}
                                                 {$rotatable_bond_name} ) {
                $rotatable_bond_name =
                    $explicit_dihedral_names->{$residue_name}
                                              {$rotatable_bond_name};
            } elsif( exists $explicit_dihedral_names->{$residue_name}
                                                      {$simplified_rotatable_bond_name} ) {
                $rotatable_bond_name =
                    $explicit_dihedral_names->{$residue_name}
                                              {$simplified_rotatable_bond_name};
            }

            $named_rotatable_bonds{$atom_id}{$rotatable_bond_name} = {
                'order' => $order{$bond_atom_ids->[1]}{$bond_atom_ids->[2]},
                'atom_ids' => $bond_atom_ids,
            };
        }
    }

    return \%named_rotatable_bonds;
}

#
# Identifies bonds that can be stretched.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $start_atom_ids - starting atom ids.
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
    my ( $bond_paths, $include_mainchain, $include_hetatoms ) = (
        $options->{'bond_paths'},
        $options->{'include_mainchain'},
        $options->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $ignore_connections = {
        'label_atom_id' => {
            'N' => { 'C' => 1 }, # Pseudo connection for heteroatoms.
        },
    };

    $bond_paths //= BondPath->new( {
        'atom_site' => $atom_site,
        'start_atom_ids' => $start_atom_ids,
        'include_hetatoms' => $include_hetatoms,
        'ignore_connections' => $ignore_connections,
    } );

    my %stretchable_bonds = ();
    my %parent_atom_ids = ();
    my %order = ();
    for my $order ( sort { $a <=> $b } keys %{ $bond_paths } ) {
        my ( $first_atom_id ) = keys %{ $bond_paths->{$order} };
        my $second_atom_id = $bond_paths->{$order}{$first_atom_id};

        $parent_atom_ids{$second_atom_id} = $first_atom_id;
        $order{$first_atom_id}{$second_atom_id} = $order;

        # Checks for mainchains and heteroatoms.
        next if ! $include_mainchain &&
            ! contains_sidechain_atoms( $parameters,
                                        $atom_site,
                                        [ $first_atom_id, $second_atom_id ] ) &&
            ! contains_hetatoms( $atom_site,
                                 [ $first_atom_id, $second_atom_id ] );

        push @{ $stretchable_bonds{$second_atom_id} },
            [ $first_atom_id, $second_atom_id ];

        # Adds bond if it is a continuation of identified bonds.
        if( exists $parent_atom_ids{$first_atom_id} &&
            exists $stretchable_bonds{$parent_atom_ids{$first_atom_id}} ) {
            unshift @{ $stretchable_bonds{$second_atom_id} },
                @{ $stretchable_bonds{$parent_atom_ids{$first_atom_id}} };
        }
    }

    # Naming the stretchable bonds.
    my %named_stretchable_bonds = ();
    for my $atom_id ( keys %stretchable_bonds ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $stretchable_bonds{$atom_id} } ) {
            my $stretchable_bond_name =
                join '-', map { $atom_site->{$_}{'label_atom_id'} }
                             @{ $bond_atom_ids };

            $named_stretchable_bonds{$atom_id}{$stretchable_bond_name} = {
                'order' => $order{$bond_atom_ids->[0]}{$bond_atom_ids->[1]},
                'atom_ids' => $bond_atom_ids,
            };
        }
    }

    return \%named_stretchable_bonds;
}

#
# Identifies bonds that can be bent.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $start_atom_ids - starting atom ids.
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
    my ( $bond_paths, $include_mainchain, $include_hetatoms ) = (
        $options->{'bond_paths'},
        $options->{'include_mainchain'},
        $options->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    $bond_paths //= BondPath->new( {
        'atom_site' => $atom_site,
        'start_atom_ids' => $start_atom_ids,
        'include_hetatoms' => $include_hetatoms,
    } );

    my %bendable_angles = ();
    my %parent_atom_ids = ();
    my %order = ();
    for my $order ( sort { $a <=> $b } keys %{ $bond_paths } ) {
        my ( $second_atom_id ) = keys %{ $bond_paths->{$order} };
        my $third_atom_id = $bond_paths->{$order}{$second_atom_id};

        $parent_atom_ids{$third_atom_id} = $second_atom_id;

        my $first_atom_id = $parent_atom_ids{$second_atom_id};

        next if ! defined $first_atom_id;

        $order{$first_atom_id}{$second_atom_id}{$third_atom_id} = $order;

        # Checks for mainchains and heteroatoms.
        next if ! $include_mainchain &&
            ! contains_sidechain_atoms( $parameters,
                                        $atom_site,
                                        [ $first_atom_id, $second_atom_id,
                                          $third_atom_id ] ) &&
            ! contains_hetatoms( $atom_site,
                                 [ $first_atom_id, $second_atom_id,
                                   $third_atom_id ] );

        push @{ $bendable_angles{$third_atom_id} },
            [ $first_atom_id, $second_atom_id, $third_atom_id ];

        # Adds bond if it is a continuation of identified bonds.
        if( exists $bendable_angles{$second_atom_id} ) {
            unshift @{ $bendable_angles{$third_atom_id} },
                @{ $bendable_angles{$second_atom_id} };
        }
    }

    my %named_bendable_angles = ();
    for my $atom_id ( keys %bendable_angles ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $bendable_angles{$atom_id} } ) {
            my $bendable_angle_name =
                join '-', map { $atom_site->{$_}{'label_atom_id'} }
                             @{ $bond_atom_ids };

            $named_bendable_angles{$atom_id}{$bendable_angle_name} = {
                'order' => $order{$bond_atom_ids->[0]}
                                 {$bond_atom_ids->[1]}
                                 {$bond_atom_ids->[2]},
                'atom_ids' => $bond_atom_ids,
            };
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
    my ( $parameters, $atom_site ) = @_;

    my $rotatable_bonds = rotatable_bonds( $parameters, $atom_site );

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

sub contains_sidechain_atoms
{
    my ( $paramters, $atom_site, $atom_ids ) = @_;
    return scalar( grep { $paramters->{'_[local]_mainchain_atom_names_table'}{$_} }
                   map  { $atom_site->{$_}{'label_atom_id'} }
                       @{ $atom_ids } ) <
           scalar( @{ $atom_ids } );
}

sub contains_hetatoms
{
    my ( $atom_site, $atom_ids ) = @_;
    return grep { $atom_site->{$_}{'group_PDB'} eq 'HETATM' } @{ $atom_ids };
}

1;
