package ConnectAtoms;

use strict;
use warnings;

BEGIN{
use Exporter qw( import );
our @EXPORT_OK = qw( append_connections
                     connect_atoms
                     connect_hetatoms
                     connect_atoms_explicitly
                     connect_two_atoms
                     is_connected
                     is_neighbour
                     is_second_neighbour
                     remove_connections
                     retains_connections );
}

use Carp qw( confess );
use Digest::MD5 qw( md5_hex );
use List::Util qw( any );

use Grid qw( identify_neighbour_cells
             grid_box );
use Measure qw( around_distance
                distance_squared );
use PDBxParser qw( filter
                   filter_new
                   split_by
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
#     $options->{'assign_hetatoms'} - if on, hetatoms are included.
# Output:
#     $is_connected - boolean: 0 (for not connected) or 1 (for connected).
#

sub is_connected
{
    my ( $parameters, $target_atom, $neighbour_atom, $options ) = @_;

    my ( $no_connection_list, $no_covalent_radii, $assign_hetatoms ) =
        ( $options->{'no_connection_list'},
          $options->{'no_covalent_radii'},
          $options->{'assign_hetatoms'} );
    $no_connection_list //= 0;
    $no_covalent_radii //= 0;
    $assign_hetatoms //= 0;

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

    if( ! exists $atom_site->{"$target_atom_id"} ) {
        confess "atom with id $target_atom_id does not exist in the atom site";
    }
    if( ! exists $atom_site->{"$target_atom_id"}{'connections'} ) {
        confess "atom with id $target_atom_id does not have 'connection' key";
    }

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

    if( ! exists $atom_site->{"$target_atom_id"} ) {
        confess "atom with id $target_atom_id does not exist in the atom site";
    }
    if( ! exists $atom_site->{"$target_atom_id"}{'connections'} ) {
        confess "atom with id $target_atom_id does not have 'connection' key";
    }

    my $is_sec_neighbour = 0;
    foreach my $i ( @{ $atom_site->{"$target_atom_id"}{'connections'} } ) {
        if( ! exists $atom_site->{$i}{'connections'} ) {
            confess "atom with id $i does not have 'connection' key";
        }

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
# Connects two atoms by including ids to 'connections' attribute in both atom
# data structure.
# Input:
#     $atom_site - $atom_site - atom site data structure (see PDBxParser.pm);
#     ${first,second}_atom_id - atom ids.
# Output:
#     none
#

sub connect_two_atoms
{
    my ( $parameters, $atom_site, $first_atom_id, $second_atom_id ) = @_;

    if( is_connected( $parameters,
                      $atom_site->{"$first_atom_id"},
                      $atom_site->{"$second_atom_id"} ) ) {
        push @{ $atom_site->{$first_atom_id}{'connections'} }, "$second_atom_id";
        push @{ $atom_site->{$second_atom_id}{'connections'} }, "$first_atom_id";
    }

    return;
}

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
#     $options->{'assign_hetatoms'} - includes hetatoms in connection
#     calculations.
# Output:
#     none - connects atoms by adding "connection" key and values to atom site
#     data structure.
#

sub connect_atoms
{
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $append_connections, $no_connection_list, $no_covalent_radii,
         $assign_hetatoms ) = (
        $options->{'append_connections'},
        $options->{'no_connection_list'},
        $options->{'no_covalent_radii'},
        $options->{'assign_hetatoms'}
    );

    $append_connections //= 0;
    $no_connection_list //= 0;
    $no_covalent_radii //= 0;
    $assign_hetatoms //= 0;

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
                                          $no_covalent_radii,
                                      'assign_hetatoms' =>
                                          $assign_hetatoms } ) ) &&
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

# NOTE: for now, hetatoms will be usually metal ions.
sub connect_hetatoms
{
    my ( $parameters, $atom_site, $options ) = @_;
    my ( $struct_conn ) = ( $options->{'struct_conn'} );

    # Connecting atoms with heteroatoms using $struct_conn.
    # NOTE: what about different models? PDBx's '_struct_conn' does not have
    # 'ptnr2_pdbx_PDB_model_num' and 'ptnr2_label_alt_id'.
    # TODO: split_by() should be more generalized. But at the moment, the
    # '_atom_site' will be used as an option.
    my $struct_conn_groups =
        split_by( { 'atom_site' => $struct_conn,
                    'attributes' => [ 'ptnr2_label_seq_id','ptnr2_label_comp_id',
                                      'ptnr2_label_asym_id' ] } );

    # TODO: do not forget to filter out non-ionic/metal atoms.
    for my $struct_conn_key ( keys %{ $struct_conn_groups } ) {
        for my $struct_conn_id ( @{ $struct_conn_groups->{$struct_conn_key} } ) {
            my ($ptnr1_label_seq_id, $ptnr1_label_comp_id, $ptnr1_label_asym_id,
                $ptnr1_label_atom_id, $ptnr1_label_alt_id,
                $ptnr2_label_seq_id, $ptnr2_label_comp_id, $ptnr2_label_asym_id,
                $ptnr2_label_atom_id, $ptnr2_label_alt_id )=
                map { $struct_conn->{$struct_conn_id}{$_} }
                    ( 'ptnr1_label_seq_id',  'ptnr1_label_comp_id',
                      'ptnr1_label_asym_id', 'ptnr1_label_atom_id',
                      'ptnr1_label_alt_id',  'ptnr2_label_seq_id',
                      'ptnr2_label_comp_id', 'ptnr2_label_asym_id',
                      'ptnr2_label_atom_id', 'ptnr2_label_alt_id' );

            next if any { ! defined $_ }
                        ( $ptnr1_label_seq_id, $ptnr1_label_comp_id,
                          $ptnr1_label_asym_id, $ptnr1_label_atom_id,
                          $ptnr2_label_seq_id, $ptnr2_label_comp_id,
                          $ptnr2_label_asym_id, $ptnr2_label_atom_id );

            # TODO: could be optimized by just storing already used hetatom
            # coordinates.
            # TODO: check if these data items specify exact atoms.
            my $hetatom_id =
                filter_new( $atom_site,
                            { 'include' =>
                              { 'label_seq_id'  => [ $ptnr2_label_seq_id ],
                                'label_comp_id' => [ $ptnr2_label_comp_id ],
                                'label_asym_id' => [ $ptnr2_label_asym_id ],
                                'label_atom_id' => [ $ptnr2_label_atom_id ],
                                'label_alt_id'  => [ defined $ptnr2_label_alt_id?
                                                     $ptnr2_label_alt_id :
                                                     '.' ] },
                               'return_data' => 'id' } );
            my $residue_atom_id =
                filter_new( $atom_site,
                            { 'include' =>
                              { 'label_seq_id'  => [ $ptnr1_label_seq_id ],
                                'label_comp_id' => [ $ptnr1_label_comp_id ],
                                'label_asym_id' => [ $ptnr1_label_asym_id ],
                                'label_atom_id' => [ $ptnr1_label_atom_id ],
                                'label_alt_id'  => [ defined $ptnr1_label_alt_id?
                                                     $ptnr1_label_alt_id :
                                                     '.' ] },
                              'return_data' => 'id' } );

            # NOTE: for now, if multiple residue atoms are present, it should be
            # skipped and standard connection with 'CA' will be performed.
            next if ! @{ $hetatom_id }      || scalar( @{ $hetatom_id } ) > 1 ||
                    ! @{ $residue_atom_id } || scalar( @{ $residue_atom_id } )>1;

            connect_atoms_explicitly( $atom_site, $hetatom_id, $residue_atom_id);
        }
    }

    # If hetatoms have no explicit connection, then they are connected to the
    # closest CAs.
    my $interaction_distance =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};
    my $interaction_atom_site =
        filter_new( $atom_site, { 'include' => { 'label_atom_id' => [ 'CA' ] }});

    my $hetatom_names =
        $parameters->{'_[local]_sidechain_hetatom_extension'};
    my $hetatom_ids =
        filter_new( $atom_site,
                    { 'include' => { 'label_comp_id' => $hetatom_names },
                          'return_data' => 'id' } );

    for my $hetatom_id ( @{ $hetatom_ids } ) {
        around_distance( $parameters,
                         $interaction_atom_site,
                         { 'id' => [ $hetatom_id ] },
                         $interaction_distance );
    }

    return;
}

# Adds connections explicitly.
# Input:
#     $atom_site - atom data structure.
#     $first_atom_id_list - first atom id list.
#     $second_atom_id_list - second atom id list.
# Output:
#     none - connects atoms by adding "connection" key and values to atom site
#     data structure.

sub connect_atoms_explicitly
{
    my ( $atom_site, $first_atom_id_list, $second_atom_id_list ) = @_;

    for my $first_atom_id ( @{ $first_atom_id_list } ) {
        for my $second_atom_id ( @{ $second_atom_id_list } ) {
            push @{ $atom_site->{$first_atom_id}{'connections'} },
                "$second_atom_id";
            push @{ $atom_site->{$second_atom_id}{'connections'} },
                "$first_atom_id";
        }
    }

    return;
}

# Removes connections explicitly.
# Input:
#     $atom_site - atom data structure.
#     $first_atom_id_list - first atom id list.
#     $second_atom_id_list - second atom id list.
# Output:
#     none - removes atom connections by deleting "connection" key and values to
#     atom site data structure.

sub remove_connections
{
    my ( $atom_site, $first_atom_id_list, $second_atom_id_list ) = @_;

    for my $first_atom_id ( @{ $first_atom_id_list } ) {
        for my $second_atom_id ( @{ $second_atom_id_list } ) {
            my ( $first_atom_id_idx ) =
                grep { $first_atom_id eq
                       $atom_site->{$second_atom_id}{'connections'}[$_]  }
                     (0..$#{ $atom_site->{$second_atom_id}{'connections'} });
            my ( $second_atom_id_idx ) =
                grep { $second_atom_id eq
                       $atom_site->{$first_atom_id}{'connections'}[$_]  }
                     (0..$#{ $atom_site->{$first_atom_id}{'connections'} });

            if( defined $second_atom_id_idx ) {
                splice @{ $atom_site->{$first_atom_id}{'connections'} },
                    $second_atom_id_idx, 1;
            }
            if( defined $first_atom_id_idx ) {
                splice @{ $atom_site->{$second_atom_id}{'connections'} },
                    $first_atom_id_idx, 1;
            }
        }
    }

    return;
}

#
# Appends one's atom ids to the others list of connections.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $appendable_site - atom site that should be appended.
# Output:
#     none
#

sub append_connections
{
    my ( $atom_site, $appendable_site ) = @_;

    for my $appendable_id ( keys %{ $appendable_site } ) {
        for my $target_atom_id ( @{ $appendable_site->{$appendable_id}
                                                      {'connections'} } ) {
            if( ! any { $appendable_id eq $_ }
                     @{ $atom_site->{$target_atom_id}{'connections'} } ) {
                push @{ $atom_site->{$target_atom_id}{'connections'} },
                    $appendable_id;
            }
        }
    }

    return;
}

sub retains_connections
{
    my ( $parameters, $atom_site_1, $atom_site_2, $options ) = @_;

    # Generates hash that translates original atom ids to atom ids from second
    # atom site. It is needed, because after rotation the atom id can change.
    my %atom_id_2_to_original = ();
    for my $atom_id_2 ( keys %{ $atom_site_2 } ) {
        my $original_atom_id_2 = original_atom_id( $atom_site_2, $atom_id_2 );
        $atom_id_2_to_original{$atom_id_2} = $original_atom_id_2;
    }

    my %atom_2_original_connections = ();
    for my $atom_id_2 ( keys %{ $atom_site_2 } ) {
        my $atom_2_connections = $atom_site_2->{$atom_id_2}{'connections'};
        if( defined $atom_2_connections ) {
            # TODO: check the correctness of the code if atom ids are missing
            # in atom site.
            $atom_2_original_connections{$atom_id_2_to_original{$atom_id_2}} =
                [ map { $atom_id_2_to_original{$_} } @{ $atom_2_connections } ];
        }
    }

    for my $atom_id_1 ( keys %{ $atom_site_1 } ) {
        return 0 if ! exists $atom_id_2_to_original{$atom_id_1};

        my $atom_1_connections = $atom_site_1->{$atom_id_1}{'connections'};

        next if ! defined $atom_1_connections;

        for my $connection_id_1 ( @{ $atom_1_connections } ) {
            return 0 if ! any { $connection_id_1 eq $_ }
                             @{ $atom_2_original_connections{$atom_id_1} };
        }
    }

    return 1;
}

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

1;
