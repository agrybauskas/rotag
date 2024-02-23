package ConnectAtoms;

use strict;
use warnings;

BEGIN{
use Exporter qw( import );
our @EXPORT_OK = qw( assign_hetatoms
                     connect_atoms
                     connect_atoms_explicitly
                     create_hetatom_struct_conn
                     is_connected
                     is_neighbour
                     is_second_neighbour );
}

use Carp qw( confess );
use Clone qw( clone );
use Digest::MD5 qw( md5_hex );
use List::Util qw( any
                   max );

use Grid qw( identify_neighbour_cells
             grid_box );
use Measure qw( around_distance
                distance_squared );
use PDBxParser qw( filter
                   filter_by_unique_residue_key
                   filter_new
                   replace_atom_site_ids
                   unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ---------------------------- Pairwise relations ----------------------------- #

#
# Checks if atoms are connected.
# Input:
#     ${target,neighbour}_atom - atom data structure (see PDBxParser.pm).
#     $options->{'no_connection_list'} - if on, connection is determined without
#     ForceField::Parameters connection list.
#     $options->{'no_covalent_radii'} - if on, connection is determined without
#     covalent radii.
# Output:
#     $is_connected - boolean: 0 (for not connected) or 1 (for connected).
#

sub is_connected
{
    my ( $parameters, $target_atom, $neighbour_atom, $options ) = @_;

    my ( $no_connection_list, $no_covalent_radii ) =
        ( $options->{'no_connection_list'},
          $options->{'no_covalent_radii'} );
    $no_connection_list //= 0;
    $no_covalent_radii //= 0;

    my $covalent_bond_comb = $parameters->{'_[local]_covalent_bond_combinations'};
    my $connectivity = $parameters->{'_[local]_connectivity'};

    # Generates all possible combinations of covalent distances.
    my $bond_length_comb =
        $covalent_bond_comb->{$target_atom->{'type_symbol'}}
                             {$neighbour_atom->{'type_symbol'}}
                             {'length'};

    my $length_error_comb =
        $covalent_bond_comb->{$target_atom->{'type_symbol'}}
                             {$neighbour_atom->{'type_symbol'}}
                             {'error'};

    if( ! defined $bond_length_comb || ! defined $length_error_comb ) {
        confess sprintf 'bond between %s and %s atoms is not characterized ' .
                        'in the current force field description',
                        $target_atom->{'type_symbol'},
                        $neighbour_atom->{'type_symbol'};
    }

    my $target_residue_name = $target_atom->{'label_comp_id'};
    my $target_atom_name = $target_atom->{'label_atom_id'};
    my $neighbour_atom_name = $neighbour_atom->{'label_atom_id'};

    # Determines their residue unique keys.
    # TODO: should be cautious about alt ids that can divide sidechain in two
    # parts.
    my $target_residue_key = unique_residue_key( $target_atom );
    my $neighbour_residue_key = unique_residue_key( $neighbour_atom );

    # Precalculates squared distance between atom pairs.
    my $distance_squared = distance_squared( $target_atom, $neighbour_atom );

    # Checks, if distance between atom pairs is in one of the combinations.
    my $bond_length;
    my $length_error;

    for my $i ( 0..$#{ $bond_length_comb } ) {
        next if $bond_length_comb->[$i][0] eq '.' ||
                $bond_length_comb->[$i][1] eq '.';

        # NOTE: connections between hetero- and other atoms have pseudo
        # connections for now.
        next if ( $target_atom->{'group_PDB'} eq 'HETATM' &&
                  $neighbour_atom->{'group_PDB'} ne 'HETATM' ) ||
            ( $target_atom->{'group_PDB'} ne 'HETATM' &&
                  $neighbour_atom->{'group_PDB'} eq 'HETATM' );

        $bond_length = $bond_length_comb->[$i][0] + $bond_length_comb->[$i][1];
        $length_error = $length_error_comb->[$i][0] + $length_error_comb->[$i][1];
        if( ( ! $no_covalent_radii ) &&
            ( $distance_squared >= ( $bond_length - $length_error ) ** 2 ) &&
            ( $distance_squared <= ( $bond_length + $length_error ) ** 2 ) ) {
            return 1;
        }
        if( ( ! $no_connection_list ) &&
            ( $target_residue_key eq $neighbour_residue_key ) &&
            ( exists $connectivity->{$target_residue_name}
                                    {$target_atom_name} &&
              any { $neighbour_atom_name eq $_  }
                 @{ $connectivity->{$target_residue_name}
                                   {$target_atom_name} } ) ) {
            return 1;
        }
    }

    return 0;
}

#
# Checks if two atoms are connected. Similar to is_connected function, but looks
# for existing "connection" keys in atom site data structure.
# Input:
#     $atom_site - atom site data structure (see PDBxParser).
#     $target_atom_id - first atom id.
#     $neighbour_atom_id - second atom id.
# Output:
#     $is_neighbour - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_neighbour
{
    my ( $atom_site, $target_atom_id, $neighbour_atom_id ) = @_;

    return 0 if ! exists $atom_site->{"$target_atom_id"} ||
        ! exists $atom_site->{"$target_atom_id"}{'connections'};

    my $is_neighbour = 0;
    foreach my $i ( @{ $atom_site->{"$target_atom_id"}{'connections'} } ) {
        if( "$neighbour_atom_id" eq "$i" ) {
            $is_neighbour = 1;
            last;
        }
    }

    return $is_neighbour;
}

#
# Checks if two atoms are separated only by one atom.
# Input:
#     $atom_site - atom site data structure.
#     $target_atom_id - id of first atom.
#     $sec_neighbour_atom_id - id of second atom.
# Output:
#     $is_sec_neighbour - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_second_neighbour
{
    my ( $atom_site, $target_atom_id, $sec_neighbour_atom_id ) = @_;

    return 0 if ! exists $atom_site->{"$target_atom_id"} ||
        ! exists $atom_site->{"$target_atom_id"}{'connections'};

    my $is_sec_neighbour = 0;
    foreach my $i ( @{ $atom_site->{"$target_atom_id"}{'connections'} } ) {
        next if ! exists $atom_site->{$i}{'connections'};

        foreach my $j ( @{ $atom_site->{$i}{'connections'} } ) {
            if( "$sec_neighbour_atom_id" eq "$j" ) {
                $is_sec_neighbour = 1;
                last;
            }
        }
        last if $is_sec_neighbour == 1;
    }

    return $is_sec_neighbour;
}

# ------------------------------ Connect atoms -------------------------------- #

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Then, all atoms' distances are compared pairwisely in 1 + 26 surrounding
# boxes. If distance correspond to appropriate length, then connection is
# made by two atoms.
# Input:
#     $atom_site - atom data structure.
#     $options->{'append_connections'} - only appends connections and not
#     recalculates the old ones.
#     $options->{'no_connection_list'} - tries to connect atoms without prior
#     knowledge about possible connections.
#     $options->{'no_covalent_radii'} - tries to connect atoms only using list
#     of connections.
# Output:
#     none - connects atoms by adding "connection" key and values to atom site
#     data structure.
#

sub connect_atoms
{
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $append_connections, $no_connection_list, $no_covalent_radii ) = (
        $options->{'append_connections'},
        $options->{'no_connection_list'},
        $options->{'no_covalent_radii'},
    );

    $append_connections //= 0;
    $no_connection_list //= 0;
    $no_covalent_radii //= 0;

    # Removes all previously described connections if certain flags are not on.
    if( ! $append_connections ) {
        for my $atom_id ( keys %{ $atom_site } ) {
            if( exists $atom_site->{$atom_id}{'connections'} ) {
                delete $atom_site->{$atom_id}{'connections'}
            }
        }
    }

    # For each cell, checks neighbouring cells. Creates box around atoms, makes
    # grid with edge length of max covalent radii of the parameter file.
    my ( $grid_box ) = grid_box( $parameters, $atom_site );
    my $neighbour_cells = identify_neighbour_cells( $grid_box );

    foreach my $cell ( keys %{ $grid_box } ) {
        foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
            foreach my $neighbour_id ( @{ $neighbour_cells->{$cell} } ) {
                if( ( is_connected( $parameters,
                                    $atom_site->{"$atom_id"},
                                    $atom_site->{"$neighbour_id"},
                                    { 'no_connection_list' =>
                                          $no_connection_list,
                                      'no_covalent_radii' =>
                                          $no_covalent_radii } ) ) &&
                    ( ( ! exists $atom_site->{$atom_id}{'connections'} ) ||
                      ( ! any { $neighbour_id eq $_ }
                             @{ $atom_site->{$atom_id}{'connections'} } ) ) ){
                    push @{ $atom_site->{$atom_id}{'connections'} },
                        "$neighbour_id";
                }
            }
        }
    }

    return;
}

# Adds connections explicitly.
# Input:
#     $atom_site - atom data structure.
#     $first_atom_id_list - first atom id list.
#     $second_atom_id_list - second atom id list.
#     $options->{'connection_type'} - connection type
#     (connections|connections_hetatom).
# Output:
#     none - connects atoms by adding "connection" key and values to atom site
#     data structure.

sub connect_atoms_explicitly
{
    my ( $atom_site, $first_atom_id_list, $second_atom_id_list, $options ) = @_;
    my ( $connection_type ) = ( $options->{'connection_type'} );

    for my $first_atom_id ( @{ $first_atom_id_list } ) {
        next if ! exists $atom_site->{$first_atom_id};

        my $first_atom_group_PDB = $atom_site->{$first_atom_id}{'group_PDB'};
        for my $second_atom_id ( @{ $second_atom_id_list } ) {
            next if ! exists $atom_site->{$second_atom_id};

            my $second_atom_group_PDB =
                $atom_site->{$second_atom_id}{'group_PDB'};
            if( ( $first_atom_group_PDB eq 'HETATM' ||
                  $second_atom_group_PDB eq 'HETATM' ) &&
                ! defined $connection_type ) {
                $connection_type = 'connections_hetatom'
            }
            $connection_type //= 'connections';
            push @{ $atom_site->{$first_atom_id}{$connection_type} },
                "$second_atom_id";
            push @{ $atom_site->{$second_atom_id}{$connection_type} },
                "$first_atom_id";
        }
    }

    return;
}

# Returns original atom id.
# Input:
#     $atom_site - atom data structure;
#     $atom_id - atom id.
# Output:
#     $original_atom_id - return original atom id.

sub original_atom_id
{
    my ( $atom_site, $atom_id ) = @_;
    my $original_atom_id = $atom_site->{$atom_id}{'original_atom_id'};
    if( defined $original_atom_id ) {
        return $original_atom_id;
    } else {
        return $atom_id;
    }
}


#
# Assigns heteroatoms to specific residues according to "_struct_conn" -- either
# creating the new one or assigning to the existing one. Duplicated heteroatoms
# are to be expected.
# Input:
#     $parameters - general parameters (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $struct_conn - reads 'struc_conn' and assings connections appropriately.
# Output:
#     atom site with assigned heteroatoms.
#

sub assign_hetatoms
{
    my ( $parameters, $atom_site, $struct_conn, $options ) = @_;

    $struct_conn //= create_hetatom_struct_conn( $parameters, $atom_site );

    return if ! %{ $struct_conn };

    my %origin_atom_site = %{ clone $atom_site };
    my %track_renamed_atom_ids = ();
    my $last_atom_id = max( keys %{ $atom_site } ) + 1;
    for my $struct_conn_id ( sort keys %{ $struct_conn } ) {
        my %related_atom_selection_1 = (
            $struct_conn->{$struct_conn_id}{'ptnr2_auth_seq_id'} eq '.' ?
            ( 'auth_seq_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr2_auth_seq_id'} ],
              'auth_asym_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr2_auth_asym_id'} ] ) :
            ( 'label_seq_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr2_label_seq_id'} ],
              'label_asym_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr2_label_asym_id'} ] )
        );

        # Related atom site is used as the connected molecule will inherit
        # some data items.
        my $related_atom_site_1 =
            filter_new( \%origin_atom_site,
                        { 'include' => \%related_atom_selection_1 } );
        my $connected_atom_site_1 = filter_new(
            $related_atom_site_1,
            { 'include' =>
              { 'label_atom_id' => [
                    $struct_conn->{$struct_conn_id}{'ptnr2_label_atom_id'} ] } }
        );

        next if ! %{ $connected_atom_site_1 };

        my %related_atom_selection_2 = (
            $struct_conn->{$struct_conn_id}{'ptnr1_auth_seq_id'} eq '.' ?
            ( 'auth_seq_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr1_auth_seq_id'} ],
              'auth_asym_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr1_auth_asym_id'} ] ) :
            ( 'label_seq_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr1_label_seq_id'} ],
              'label_asym_id' => [
                  $struct_conn->{$struct_conn_id}{'ptnr1_label_asym_id'} ] ),
            ( 'label_atom_id' => [
                $struct_conn->{$struct_conn_id}{'ptnr1_label_atom_id'} ] )
        );

        my $connected_atom_site_2 = filter_new(
            \%origin_atom_site,
            { 'include' =>
              { $related_atom_selection_2,
                'label_atom_id' => [
                    $struct_conn->{$struct_conn_id}{'ptnr1_label_atom_id'} ] } }
        );

        next if ! %{ $connected_atom_site_2 };

        my ( $connected_atom_id_1 ) = keys %{ $connected_atom_site_1 };
        my ( $connected_atom_id_2 ) = keys %{ $connected_atom_site_2 };

        # # Especially important when there are heteroatoms connected to another
        # # heteroatoms
        # if( exists $track_renamed_atom_ids{$connected_atom_id_1} ) {
        #     $connected_atom_id_1 =
        #         $track_renamed_atom_ids{$connected_atom_id_1};
        # }
        # if( exists $track_renamed_atom_ids{$connected_atom_id_2} ) {
        #     $connected_atom_id_2 =
        #         $track_renamed_atom_ids{$connected_atom_id_2};
        # }

        for my $related_atom_id_1 ( sort keys %{ $related_atom_site_1 } ) {
            replace_atom_site_ids( $atom_site,
                                   [ { 'from' => $related_atom_id_1,
                                       'to' => $last_atom_id } ],
                                   $options );

            $track_renamed_atom_ids{$related_atom_id_1} = $last_atom_id;

            if( $related_atom_id_1 eq $connected_atom_id_1 ) {
                connect_atoms_explicitly( $atom_site,
                                          [ $last_atom_id ],
                                          [ $connected_atom_id_2 ] );
            }

            for my $attribute ( 'label_seq_id', 'label_asym_id', 'label_alt_id',
                                'pdbx_PDB_model_num'  ) {
                $atom_site->{$last_atom_id}{$attribute} =
                    $origin_atom_site{$connected_atom_id_2}{$attribute};
            }

            $last_atom_id++;
        }
    }

    return;
}

#
# Creates "_struct_conn" data structure that shows explicit interactions between
# heteroatoms  specific residues connecting them to N, O, P, S atoms -- either
# creating the new one or assigning to the existing one.
# Input:
#     $parameters - general parameters (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_ids - atom ids.
# Output:
#     "_struct_conn" data structure.
#

sub create_hetatom_struct_conn
{
    my ( $parameters, $atom_site ) = @_;

    my $hetatom_site =
        filter_new( $atom_site,
                    { 'include' => { 'group_PDB' => [ 'HETATM' ] } } );
    my $interaction_atom_site =
        filter_new( $atom_site,
                    { 'include' =>
                      { 'type_symbol' => [ 'H', 'N', 'O', 'S' ] } } );

    my %struct_conn = ();
    my $metalc_counter = 1;
    my $hydrog_counter = 1;
    for my $hetatom_id ( sort keys %{ $hetatom_site } ) {
        my $residue_name = $hetatom_site->{$hetatom_id}{'type_symbol'};
        my $struct_conn_type = $residue_name eq 'HOH' ? 'hydrog' : 'metalc';
        my $residue_interaction_site;
        my $interaction_distance;
        if( $struct_conn_type eq 'hydrog' ) {
            $interaction_distance =
                $parameters->{'_[local]_constants'}{'hydrog_length'}
                             {$hetatom_site->{$hetatom_id}{'type_symbol'}};
            # NOTE: study which atoms should be included -- are we looking
            # at hydrogen bonding or induced non-covalent bonding.
            $residue_interaction_site = $interaction_atom_site;
        } else {
            $interaction_distance =
                $parameters->{'_[local]_constants'}{'metalc_length'}
                             {$hetatom_site->{$hetatom_id}{'type_symbol'}};
            $residue_interaction_site =
                { map { ( $_ => $interaction_atom_site->{$_} ) }
                  grep { $interaction_atom_site->{$_}{'type_symbol'} ne 'H'  }
                  keys %{ $interaction_atom_site } };
        }

        next if ! defined $interaction_distance;

        my $around_site =
            around_distance( $parameters,
                             { $hetatom_id => $hetatom_site->{$hetatom_id},
                               %{ $interaction_atom_site } },
                             { 'id' => [ $hetatom_id ] },
                             $interaction_distance );

        next if ! %{ $around_site };

        for my $around_atom_id ( sort { $a <=> $b } keys %{ $around_site } ) {
            next if $atom_site->{$around_atom_id}{'type_symbol'} eq 'H' &&
                ! any { $atom_site->{$_}{'type_symbol'} eq 'N' ||
                        $atom_site->{$_}{'type_symbol'} eq 'O' ||
                        $atom_site->{$_}{'type_symbol'} eq 'S' }
                     @{ $atom_site->{$around_atom_id}{'connections'} };

            my $struct_conn_id =
                $struct_conn_type .
                ( $struct_conn_type eq 'metalc' ?
                  $metalc_counter : $hydrog_counter );

            $struct_conn{$struct_conn_id} = {
                'id' => $struct_conn_id,
                'conn_type_id' => $struct_conn_type,
                'ptnr1_label_atom_id' =>
                    $atom_site->{$around_atom_id}{'label_atom_id'},
                'ptnr1_label_asym_id' =>
                    $atom_site->{$around_atom_id}{'label_asym_id'},
                'ptnr1_label_comp_id' =>
                    $atom_site->{$around_atom_id}{'label_comp_id'},
                'ptnr1_label_seq_id' =>
                    $atom_site->{$around_atom_id}{'label_seq_id'},
                'ptnr2_label_atom_id' =>
                    $atom_site->{$hetatom_id}{'label_atom_id'},
                'ptnr2_label_asym_id' =>
                    $atom_site->{$hetatom_id}{'label_asym_id'},
                'ptnr2_label_comp_id' =>
                    $atom_site->{$hetatom_id}{'label_comp_id'},
                'ptnr2_label_seq_id' =>
                    $atom_site->{$hetatom_id}{'label_seq_id'},
            };

            $metalc_counter++ if $struct_conn_type eq 'metalc';
            $hydrog_counter++ if $struct_conn_type eq 'hydrog';
        }
    }

    return \%struct_conn;
}

1;
