package BondParameters;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( bendable_angles
                     collect_bond_angles
                     collect_bond_lengths
                     collect_dihedral_angles
                     combine_bond_parameters
                     restructure_by_atom_ids
                     rotatable_bonds
                     stretchable_bonds );

use Carp;
use List::Util qw( any
                   none );

use BondPath;
use AtomProperties qw( contains_hetatoms
                       contains_sidechain_atoms );
use Measure qw( bond_angle
                bond_length
                dihedral_angle );
use PDBxParser qw( expand
                   filter_new
                   split_by );

#
# Identifies bonds that can be rotated by torsional angle.
# Input:
#     $parameters - parameter data structure (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $options->{'start_atom_ids'} - starting atom ids;
#     $options->{'include_mainchain'} - flag that includes main-chain atoms;
#     $options->{'include_hetatoms'} - flag that includes heteroatoms.
# Output:
#     adds data structure that describes rotatable bonds and the constituent atom
#     ids of the bond.
#

sub rotatable_bonds
{
    my ( $parameters, $atom_site, $options ) = @_;
    my ( $start_atom_ids, $include_mainchain, $include_hetatoms ) = (
        $options->{'start_atom_ids'},
        $options->{'include_mainchain'},
        $options->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $explicit_dihedral_names = $parameters->{'_[local]_dihedral_angle_name'};
    my $ignore_connections = {
        'label_atom_id' => {
            'N' => { 'CD' => 1 },
        },
    };

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %rotatable_bonds = ();
    my %bond_order = ();
    my $bond_order_idx = 1;

    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $residue_site =
            filter_new( $atom_site,
                        { 'include' =>
                          { 'id' => $residue_groups->{$residue_unique_key} } } );

        if( $include_mainchain ) {
            my @expanded_atom_ids = @{ expand( $residue_site, $atom_site, 1 ) };
            my ( $residue_id, $chain_id, $pdbx_model_id, $alt_id ) =
                split ',', $residue_unique_key;
            $start_atom_ids =
                filter_new( $residue_site,
                            { 'include' => { 'id' => \@expanded_atom_ids,
                                             'label_atom_id' => [ 'C' ],
                                             'label_asym_id' => [ $chain_id ],
                                             'pdbx_PDB_model_num' => [ $pdbx_model_id ],
                                             'label_alt_id' => [ $alt_id ] },
                              'exclude' => { 'label_seq_id' => [ $residue_id ] },
                              'return_data' => 'id' } );
            $start_atom_ids = @{ $start_atom_ids } ? $start_atom_ids : undef;
        }

        my $bond_paths = BondPath->new( {
            'atom_site' => $residue_site,
            'start_atom_ids' => $start_atom_ids,
            'include_hetatoms' => $include_hetatoms,
            'ignore_connections' => $ignore_connections,
        } );

        my %rotatable_bonds_cache = ();
        for my $fourth_atom_id ( @{ $bond_paths->get_atom_order } ) {
            my $third_atom_id = $bond_paths->get_atom_id_to( $fourth_atom_id );

            next if ! defined $third_atom_id;

            my $second_atom_id = $bond_paths->get_atom_id_to( $third_atom_id );

            next if ! defined $second_atom_id;

            my $first_atom_id = $bond_paths->get_atom_id_to( $second_atom_id );

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

            # In order not to get unnecessary bond parameters, there should be
            # no bond parameters that contains pseudo connections, but do not
            # contain heteroatoms.
            next if ! contains_hetatoms( $atom_site,
                                         [ $first_atom_id, $second_atom_id,
                                           $third_atom_id, $fourth_atom_id ] ) &&
                any { $bond_paths->get_connection_type( $_->[0], $_->[1] ) eq 'connections_hetatom' }
                     ( [ $first_atom_id, $second_atom_id ],
                       [ $second_atom_id, $third_atom_id ],
                       [ $third_atom_id, $fourth_atom_id ] );

            # Check on hybridization.
            if( ! exists $atom_site->{$second_atom_id}{'hybridization'} ) {
                confess "atom with id $second_atom_id lacks information about " .
                    "hybridization";
            }
            if( ! exists $atom_site->{$third_atom_id}{'hybridization'} ){
                confess "atom with id $third_atom_id lacks information " .
                    "about hybridization";
            }

            # If one of the bond atoms are sp3, it is rotatable.
            if( $atom_site->{$second_atom_id}{'hybridization'} eq 'sp3' ||
                $atom_site->{$third_atom_id}{'hybridization'} eq 'sp3' ||
                ( $include_hetatoms &&
                  $atom_site->{$fourth_atom_id}{'group_PDB'} eq 'HETATM' &&
                  $atom_site->{$fourth_atom_id}{'hybridization'} eq '.' ) ) {
                if( exists $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} ) {
                    push @{ $rotatable_bonds{$fourth_atom_id} },
                        $rotatable_bonds_cache{$second_atom_id}{$third_atom_id};
                } else {
                    push @{ $rotatable_bonds{$fourth_atom_id} },
                        [ $first_atom_id, $second_atom_id, $third_atom_id,
                          $fourth_atom_id ];
                    $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} =
                        [ $first_atom_id, $second_atom_id, $third_atom_id,
                          $fourth_atom_id ];
                }
            }

            # If bond atoms are sp2/sp (do not rotate) or just is a continuation
            # of the bond chain, inherits its previous atom's rotatable bonds.
            if( exists $rotatable_bonds{$third_atom_id} ) {
                unshift @{ $rotatable_bonds{$fourth_atom_id} },
                    grep { scalar( grep { defined $residue_site->{$_} } @{ $_ } ) == 4 }
                        @{ $rotatable_bonds{$third_atom_id} };
            }

            $bond_order{$second_atom_id}{$third_atom_id} = $bond_order_idx;
            $bond_order_idx++;
        }
    }

    # Naming the rotatable bonds and calculating their values.
    my %dihedral_angles_cache = ();
    for my $atom_id ( keys %rotatable_bonds ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $rotatable_bonds{$atom_id} } ) {
            my @rotatable_bond_name_keys =
                ( join( '-', map { $atom_site->{$_}{'label_atom_id'} }
                                @{ $bond_atom_ids } ),
                  join( '-', ( '.',
                               $atom_site->{$bond_atom_ids->[1]}{'label_atom_id'},
                               $atom_site->{$bond_atom_ids->[2]}{'label_atom_id'},
                               '.' ) ) );
            my ( $rotatable_bond_name ) =
                grep { defined $_ }
                     ( ( map { $explicit_dihedral_names->{$residue_name}{$_} }
                             @rotatable_bond_name_keys ),
                       ( map { $explicit_dihedral_names->{'.'}{$_} }
                             @rotatable_bond_name_keys ),
                       $rotatable_bond_name_keys[0] );

            # Calculates dihedral angle if it is not already calculated.
            my $dihedral_angle_key = join ",", @{ $bond_atom_ids };
            if( ! defined $dihedral_angles_cache{$dihedral_angle_key} ) {
                $dihedral_angles_cache{$dihedral_angle_key} = dihedral_angle(
                    [ map { [ $atom_site->{$_}{'Cartn_x'},
                              $atom_site->{$_}{'Cartn_y'},
                              $atom_site->{$_}{'Cartn_z'} ] }
                         @{ $bond_atom_ids } ]
                    );
            }

            $atom_site->{$atom_id}{'rotatable_bonds'}{$rotatable_bond_name} = {
                'order' => $bond_order{$bond_atom_ids->[1]}{$bond_atom_ids->[2]},
                'atom_ids' => $bond_atom_ids,
                'value' => $dihedral_angles_cache{$dihedral_angle_key}
            };
        }
    }

    return;
}

#
# Identifies bonds that can be stretched.
# Input:
#     $parameters - parameter data structure (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $options->{'start_atom_ids'} - starting atom ids;
#     $options->{'include_mainchain'} - flag that includes main-chain atoms;
#     $options->{'include_hetatoms'} - flag that includes heteroatoms.
# Output:
#     adds data structure that describes stretchable bonds and the constituent
#     atom ids of the bond.
#

sub stretchable_bonds
{
    my ( $parameters, $atom_site, $options  ) = @_;
    my ( $start_atom_ids, $include_mainchain, $include_hetatoms ) = (
        $options->{'start_atom_ids'},
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

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %stretchable_bonds = ();
    my %bond_order = ();
    my $bond_order_idx = 1;

    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $residue_site =
            filter_new( $atom_site,
                        { 'include' =>
                          { 'id' => $residue_groups->{$residue_unique_key} } } );

        if( $include_mainchain ) {
            my @expanded_atom_ids = @{ expand( $residue_site, $atom_site, 1 ) };
            my ( $residue_id, $chain_id, $pdbx_model_id, $alt_id ) =
                split ',', $residue_unique_key;
            $start_atom_ids =
                filter_new( $residue_site,
                            { 'include' => { 'id' => \@expanded_atom_ids,
                                             'label_atom_id' => [ 'C' ],
                                             'label_asym_id' => [ $chain_id ],
                                             'pdbx_PDB_model_num' => [ $pdbx_model_id ],
                                             'label_alt_id' => [ $alt_id ] },
                              'exclude' => { 'label_seq_id' => [ $residue_id ] },
                              'return_data' => 'id' } );
            $start_atom_ids = @{ $start_atom_ids } ? $start_atom_ids : undef;
        }

        my $bond_paths = BondPath->new( {
            'atom_site' => $residue_site,
            'start_atom_ids' => $start_atom_ids,
            'include_hetatoms' => $include_hetatoms,
            'ignore_connections' => $ignore_connections,
        } );

        for my $second_atom_id ( @{ $bond_paths->get_atom_order } ) {
            my $first_atom_id = $bond_paths->get_atom_id_to( $second_atom_id );

            next if ! defined $first_atom_id;

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
            if( exists $stretchable_bonds{$first_atom_id} ) {
                unshift @{ $stretchable_bonds{$second_atom_id} },
                    grep { scalar( grep { defined $residue_site->{$_} } @{ $_ } ) == 2 }
                        @{ $stretchable_bonds{$first_atom_id} };
            }

            $bond_order{$first_atom_id}{$second_atom_id} = $bond_order_idx;
            $bond_order_idx++;
        }
    }

    # Naming the stretchable bonds and calculating their values.
    my %bond_lengths_cache = ();
    for my $atom_id ( keys %stretchable_bonds ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $stretchable_bonds{$atom_id} } ) {
            my $stretchable_bond_name =
                join '-', map { $atom_site->{$_}{'label_atom_id'} }
                             @{ $bond_atom_ids };

            # Calculates bond length if it is not already calculated.
            my $bond_length_key = join ",", @{ $bond_atom_ids };
            if( ! defined $bond_lengths_cache{$bond_length_key} ) {
                $bond_lengths_cache{$bond_length_key} = bond_length(
                    [ map { [ $atom_site->{$_}{'Cartn_x'},
                              $atom_site->{$_}{'Cartn_y'},
                              $atom_site->{$_}{'Cartn_z'} ] }
                         @{ $bond_atom_ids } ]
                );
            }

            $atom_site->{$atom_id}{'stretchable_bonds'}{$stretchable_bond_name} = {
                'order' => $bond_order{$bond_atom_ids->[0]}{$bond_atom_ids->[1]},
                'atom_ids' => $bond_atom_ids,
                'value' => $bond_lengths_cache{$bond_length_key}
            };
        }
    }

    return;
}

#
# Identifies bonds that can be bended.
# Input:
#     $parameters - parameter data structure (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $options->{'start_atom_ids'} - starting atom ids;
#     $options->{'include_mainchain'} - flag that includes main-chain atoms;
#     $options->{'include_hetatoms'} - flag that includes heteroatoms.
# Output:
#     adds data structure that describes bendable bonds and the constituent atom
#     ids of the bond.
#

sub bendable_angles
{
    my ( $parameters, $atom_site, $options ) = @_;
    my ( $start_atom_ids, $include_mainchain, $include_hetatoms ) = (
        $options->{'start_atom_ids'},
        $options->{'include_mainchain'},
        $options->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %bendable_angles = ();
    my %bond_order = ();
    my $bond_order_idx = 1;

    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $residue_site =
            filter_new( $atom_site,
                        { 'include' =>
                          { 'id' => $residue_groups->{$residue_unique_key} } } );

        if( $include_mainchain ) {
            my @expanded_atom_ids = @{ expand( $residue_site, $atom_site, 1 ) };
            my ( $residue_id, $chain_id, $pdbx_model_id, $alt_id ) =
                split ',', $residue_unique_key;
            $start_atom_ids =
                filter_new( $residue_site,
                            { 'include' => { 'id' => \@expanded_atom_ids,
                                             'label_atom_id' => [ 'C' ],
                                             'label_asym_id' => [ $chain_id ],
                                             'pdbx_PDB_model_num' => [ $pdbx_model_id ],
                                             'label_alt_id' => [ $alt_id ] },
                              'exclude' => { 'label_seq_id' => [ $residue_id ] },
                              'return_data' => 'id' } );
            $start_atom_ids = @{ $start_atom_ids } ? $start_atom_ids : undef;
        }

        my $bond_paths = BondPath->new( {
            'atom_site' => $residue_site,
            'start_atom_ids' => $start_atom_ids,
            'include_hetatoms' => $include_hetatoms,
        } );

        for my $third_atom_id ( @{ $bond_paths->get_atom_order } ) {
            my $second_atom_id = $bond_paths->get_atom_id_to( $third_atom_id );

            next if ! defined $second_atom_id;

            my $first_atom_id = $bond_paths->get_atom_id_to( $second_atom_id );

            next if ! defined $first_atom_id;

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
                    grep { scalar( grep { defined $residue_site->{$_} } @{ $_ } ) == 3 }
                        @{ $bendable_angles{$second_atom_id} };
            }

            $bond_order{$first_atom_id}{$second_atom_id}{$third_atom_id} = $bond_order_idx;
            $bond_order_idx++;
        }
    }

    # Naming the bendable angles and calculating their values.
    my %bond_angles_cache = ();
    for my $atom_id ( keys %bendable_angles ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $bendable_angles{$atom_id} } ) {
            my $bendable_angle_name =
                join '-', map { $atom_site->{$_}{'label_atom_id'} }
                             @{ $bond_atom_ids };

            # Calculates bond angle if it is not already calculated.
            my $bond_angle_key = join ",", @{ $bond_atom_ids };
            if( ! defined $bond_angles_cache{$bond_angle_key} ) {
                $bond_angles_cache{$bond_angle_key} = bond_angle(
                    [ map { [ $atom_site->{$_}{'Cartn_x'},
                              $atom_site->{$_}{'Cartn_y'},
                              $atom_site->{$_}{'Cartn_z'} ] }
                         @{ $bond_atom_ids } ]
                );
            }

            $atom_site->{$atom_id}{'bendable_angles'}{$bendable_angle_name} = {
                'order' => $bond_order{$bond_atom_ids->[0]}
                                      {$bond_atom_ids->[1]}
                                      {$bond_atom_ids->[2]},
                'atom_ids' => $bond_atom_ids,
                'value' => $bond_angles_cache{$bond_angle_key}
            };
        }
    }

    return;
}

#
# Collects dihedral angles for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
#     'rotatable_bonds' value must be present that is generated with
#     rotatable_bonds().
# Output:
#     %dihedral_angles - returns non-redundant dihedral angles.

sub collect_dihedral_angles
{
    my ( $atom_site, $options ) = @_;
    my ( $ignore_dihedral_angles ) = ( $options->{'ignore_dihedral_angles'} );

    $ignore_dihedral_angles //= {
        'N-CA-C-O' => 1,
    };

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %dihedral_angles = ();
    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $unique_rotatable_bonds = unique_bond_parameters(
            { map { ( $_ => $atom_site->{$_}{'rotatable_bonds'} ) }
              grep { defined $atom_site->{$_}{'rotatable_bonds'} }
                  @{ $residue_groups->{$residue_unique_key} } }
        );
        $dihedral_angles{$residue_unique_key} = {
            map { ( $_ => $unique_rotatable_bonds->{$_} ) }
            grep { ! ( defined $ignore_dihedral_angles->{$_} &&
                       $ignore_dihedral_angles->{$_} ) }
            keys %{ $unique_rotatable_bonds }
        };
    }

    return \%dihedral_angles;
}

#
# Collects bond lengths for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
#     'stretchable_bonds' value must be present that is generated with
#     stretchable_bonds().
# Output:
#     %bond_lengths - returns non-redundant bond lengths.

sub collect_bond_lengths
{
    my ( $atom_site ) = @_;

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %bond_lengths = ();
    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $unique_stretchable_bonds = unique_bond_parameters(
            { map { ( $_ => $atom_site->{$_}{'stretchable_bonds'} ) }
              grep { defined $atom_site->{$_}{'stretchable_bonds'} }
                  @{ $residue_groups->{$residue_unique_key} } }
        );
        $bond_lengths{$residue_unique_key} = $unique_stretchable_bonds;
    }

    return \%bond_lengths;
}

#
# Collects bond angles for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
#     'bendable_angles' value must be present that is generated with
#     bendable_angles().
# Output:
#     %bond_angles - returns non-redundant bond angles.

sub collect_bond_angles
{
    my ( $atom_site ) = @_;

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    my %bond_angles = ();
    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $unique_bendable_angles = unique_bond_parameters(
            { map { ( $_ => $atom_site->{$_}{'bendable_angles'} ) }
              grep { defined $atom_site->{$_}{'bendable_angles'} }
                  @{ $residue_groups->{$residue_unique_key} } }
        );
        $bond_angles{$residue_unique_key} = $unique_bendable_angles;
    }

    return \%bond_angles;
}

#
# Reduces bond parameters to non-redundant ones.
# Input:
#     $bond_parameters - bond parameter data structure produced by
#     rotatable_bonds(), stretchable_bonds() or bendable_angles().
# Output:
#     %unique_bond_parameters - unique bond parameter values.

sub unique_bond_parameters
{
    my ( $bond_parameters ) = @_;
    my %unique_bond_parameters;
    for my $atom_id ( sort keys %{ $bond_parameters } ) {
        for my $parameter_name ( keys %{ $bond_parameters->{"$atom_id"} } ) {
            if( ! exists $unique_bond_parameters{"$parameter_name"} ) {
                $unique_bond_parameters{"$parameter_name"} =
                    $bond_parameters->{"$atom_id"}{"$parameter_name"};
            } else {
                # Sometimes there might be bond parameter namespace clashes so,
                # the bond parameters with highest order are used.
                $unique_bond_parameters{"$parameter_name"} =
                    $bond_parameters->{"$atom_id"}{"$parameter_name"}{'order'} >
                    $unique_bond_parameters{"$parameter_name"}{'order'} ?
                    $bond_parameters->{"$atom_id"}{"$parameter_name"} :
                    $unique_bond_parameters{"$parameter_name"};
            }
        }
    }
    return \%unique_bond_parameters;
}

#
# Restructures bond parameters by atom ids.
# Input:
#     $bond_parameters - bond parameter data structure produced by
#     rotatable_bonds(), stretchable_bonds() or bendable_angles().
# Output:
#     %restructured_bond_parameters - restructure bond parameter values.

sub restructure_by_atom_ids
{
    my ( $bond_parameters ) = @_;
    my ( $unique_residue_key ) = keys %{ $bond_parameters };
    my $restructured = {};
    for my $parameter_name ( keys %{ $bond_parameters->{$unique_residue_key} } ){
        my $parameter_data =
            { %{ $bond_parameters->{$unique_residue_key}{$parameter_name} },
              'name' => $parameter_name };
        my $atom_ids = $parameter_data->{'atom_ids'};
        my $last_hash_addr = \$restructured;
        my $last_atom_id;
        for my $atom_id ( @{ $atom_ids } ) {
            if( ! defined ${ $last_hash_addr }->{$atom_id} ) {
                ${ $last_hash_addr }->{$atom_id} = {};
            }
            $last_hash_addr = \${ $last_hash_addr }->{$atom_id};
            $last_atom_id = $atom_id;
        }
        ${ $last_hash_addr } = $parameter_data;
    }
    return $restructured;
}

sub combine_bond_parameters
{
    my ( $bond_parameters ) = @_;
    my %combined_bond_parameters = ();
    for my $bond_parameter ( @{ $bond_parameters } ) {
        for my $parameter_name ( keys %{ $bond_parameter } ) {
            next if exists $combined_bond_parameters{$parameter_name};

            $combined_bond_parameters{$parameter_name} =
                $bond_parameter->{$parameter_name};
        }
    }
    return \%combined_bond_parameters;
}

1;
