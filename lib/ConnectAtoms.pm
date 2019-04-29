package ConnectAtoms;

use strict;
use warnings;

BEGIN{
use Exporter qw( import );
our @EXPORT_OK = qw( append_connections
                     connect_atoms
                     connection_digest
                     connect_two_atoms
                     is_connected
                     is_neighbour
                     is_second_neighbour
                     retains_connections );
}

use Carp qw( confess );
use Digest::MD5 qw( md5_hex );
use List::Util qw( any );

use Grid qw( identify_neighbour_cells
             grid_box );
use Measure qw( distance_squared );
use PDBxParser qw( filter
                   unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ---------------------------- Pairwise relations ----------------------------- #

#
# Checks if atoms are connected.
# Input:
#     ${target,neighbour}_atom - atom data structure (see PDBxParser.pm).
#     $options->{'only_covalent_radii'} - if on, connection is determined only
#     with covalent radii and ForceField::Parameters is ignored.
# Output:
#     $is_connected - boolean: 0 (for not connected) or 1 (for connected).
#

sub is_connected
{
    my ( $parameters, $target_atom, $neighbour_atom, $options ) = @_;

    my ( $only_covalent_radii ) = ( $options->{'only_covalent_radii'}, );
    $only_covalent_radii //= 0;

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
    my $is_connected;

    for my $i ( 0..$#{ $bond_length_comb } ) {
        $bond_length = $bond_length_comb->[$i][0] + $bond_length_comb->[$i][1];
        $length_error = $length_error_comb->[$i][0] + $length_error_comb->[$i][1];
        if( ( $distance_squared >= ( $bond_length - $length_error ) ** 2 ) &&
            ( $distance_squared <= ( $bond_length + $length_error ) ** 2 ) ) {
            $is_connected = 1;
            last;
        } elsif( ( ! $only_covalent_radii ) &&
                 ( $target_residue_key eq $neighbour_residue_key ) &&
                 ( exists $connectivity->{$target_residue_name}
                                         {$target_atom_name} &&
                   any { $neighbour_atom_name eq $_  }
                      @{ $connectivity->{$target_residue_name}
                                        {$target_atom_name} } ) ) {
            $is_connected = 1;
        } else {
            $is_connected = 0;
        }
    }

    return $is_connected;
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
# Output:
#     none - connects atoms by adding "connection" key and values to atom site
#     data structure.
#

sub connect_atoms
{
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $append_connections, $only_covalent_radii ) = (
        $options->{'append_connections'},
        $options->{'only_covalent_radii'},
    );

    $append_connections //= 0;
    $only_covalent_radii //= 0;

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
                                    { 'only_covalent_radii' =>
                                          $only_covalent_radii } ) ) &&
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

#
# Generates hash heys or digests for given atom site with connections inside.
# Input:
#     $parameter - parameter object;
#     $atom_site - atom site data structure.
# Output:
#     $connection_digest - string that summarizes the connections of given
#     atom site and returns md5 hex.
#

sub connection_digest
{
    my ( $parameters, $atom_site, $options ) = @_;

    my ( $show_string ) = $options->{'show_string'};
    $show_string //= 0;

    # Sorts atom ids by the atom ids and original atom ids if they are pseudo
    # atoms.
    # TODO: try to find a way to check even if there are multiple pseudo atoms
    # with the same original atom ids.
    my @sorted_atom_ids =
        sort  { original_atom_id( $atom_site, $a ) <=>
                original_atom_id( $atom_site, $b ) }
        keys %{ $atom_site };

    my @digest_list = ();
    for my $atom_id ( @sorted_atom_ids ) {
        my $connections = $atom_site->{$atom_id}{'connections'};
        my @atom_ids_only_original = ();

        my $original_atom_id = $atom_site->{$atom_id}{'original_atom_id'};
        if( defined $original_atom_id  ) {
            push @digest_list, $original_atom_id;
        } else {
            push @digest_list, $atom_id;
        }

        for my $connection_id ( sort @{ $connections } ) {
            my $connection_original_id =
                $atom_site->{$connection_id}{'original_atom_id'};
            if( defined $connection_original_id ) {
                push @atom_ids_only_original, $connection_original_id;
            } else {
                push @atom_ids_only_original, $connection_id;
            }
        }

        push @digest_list, sort { $a <=> $b } @atom_ids_only_original;
    }

    my $connection_digest = join ',', @digest_list;
    if( ! $show_string ) {
        $connection_digest = md5_hex( join ',', @digest_list );
    }

    return $connection_digest;
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

sub retains_connections
{

}

1;
