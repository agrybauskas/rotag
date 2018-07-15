package ConnectAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( around_distance
		     connect_atoms
		     distance
		     distance_squared
		     is_connected
		     is_neighbour
		     is_second_neighbour );

use List::Util qw( any
		   max
		   min );

use AtomProperties qw( %ATOMS );
use Combinatorics qw( permutation );
use Grid qw( create_box
             identify_neighbour_cells
             grid_box );
use PDBxParser qw( filter );

# ------------------------------ Connect atoms ------------------------------- #

sub is_connected
{
    my ( $target_atom, $neighbour_atom ) = @_;

    # Generates all possible combinations of covalent distances.
    my $bond_length_comb =
	permutation( 2,
		     [],
		     [ $ATOMS{$target_atom->{'type_symbol'}}
			     {'covalent_radius'}
			     {'length'},
		       $ATOMS{$neighbour_atom->{'type_symbol'}}
			     {'covalent_radius'}
			     {'length'} ],
		     [] );

    my $length_error_comb =
	permutation( 2,
		     [],
		     [ $ATOMS{$target_atom->{'type_symbol'}}
			     {'covalent_radius'}
			     {'error'},
		       $ATOMS{$neighbour_atom->{'type_symbol'}}
			     {'covalent_radius'}
			     {'error'} ],
		     [] );

    # Precalculates distance between atom pairs.
    my $distance =
	( $neighbour_atom->{'Cartn_x'} - $target_atom->{'Cartn_x'} ) ** 2
      + ( $neighbour_atom->{'Cartn_y'} - $target_atom->{'Cartn_y'} ) ** 2
      + ( $neighbour_atom->{'Cartn_z'} - $target_atom->{'Cartn_z'} ) ** 2;

    # Checks, if distance between atom pairs is in one of the combinations.
    my $bond_length;
    my $length_error;
    my $is_connected;

    for( my $i = 0; $i < scalar( @{ $bond_length_comb } ); $i++ ) {
	$bond_length =
	    $bond_length_comb->[$i][0] + $bond_length_comb->[$i][1];
	$length_error =
	    $length_error_comb->[$i][0] + $length_error_comb->[$i][1];
	if( ( $distance >= ( $bond_length - $length_error ) ** 2 )
	 && ( $distance <= ( $bond_length + $length_error ) ** 2 ) ) {
	    $is_connected = 1;
	    last;
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
#     $target_id - first atom id.
#     $neighbour_id - second atom id.
# Output:
#     $is_neighbour - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_neighbour
{
    my ( $atom_site, $target_atom_id, $neighbour_id ) = @_;

    my $is_neighbour = 0;
    foreach my $i ( @{ $atom_site->{"$target_atom_id"}{'connections'} } ) {
	if( "$neighbour_id" eq "$i" ) {
	    $is_neighbour = 1;
	    last;
	}
    }

    return $is_neighbour;
}

#
# Checks if two atoms are separated by one atom.
# Input:
#     $atom_site - atom data structure.
#     $target_atom_id - id of first atom.
#     $sec_neighbour_id - id of second atom.
# Output:
#     $is_sec_neighbour - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_second_neighbour
{
    my ( $atom_site, $target_atom_id, $sec_neighbour_id ) = @_;

    my $is_sec_neighbour = 0;

    foreach my $i (
	@{ $atom_site->{"$target_atom_id"}{'connections'} } ) {
    foreach my $j (
	@{ $atom_site->{$i}{'connections'} } ) {
	if( "$sec_neighbour_id" eq "$j" ) {
	    $is_sec_neighbour = 1;
	    last;
	}
    } last if $is_sec_neighbour == 1; }

    return $is_sec_neighbour;
}

sub distance_squared
{
    my ( $atom_i, $atom_j ) = @_;

    my $distance_squared =
	( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
      + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
      + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    return $distance_squared;
}

sub distance
{
    my ( $atom_i, $atom_j ) = @_;

    return sqrt( distance_squared( $atom_i, $atom_j ) );
}

sub around_distance
{
    my ( $atom_site, $atom_specifier, $distance ) = @_;

    my @atom_ids = @{ filter( { 'atom_site' => $atom_site,
				'include' => $atom_specifier,
				'data' => [ 'id' ],
				'is_list' => 1 } ) };

    # For each cell, checks neighbouring cells. Creates box around atoms, makes
    # grid with edge length of max covalent radii of the parameter file.
    my @cell_indexes;
    my ( $grid_box, $atom_cell_pos ) =
	grid_box( $atom_site, $distance * 2, \@atom_ids );
    my $neighbour_cells = identify_neighbour_cells( $grid_box, $atom_cell_pos );

    # Checks for neighbouring cells for each cell.
    my %around_atom_site;
    foreach my $cell ( keys %{ $atom_cell_pos } ) {
        foreach my $atom_id ( @{ $atom_cell_pos->{$cell} } ) {
            foreach my $neighbour_id ( @{ $neighbour_cells->{$cell} } ) {
        	if( ( ! any { $neighbour_id eq $_ } @atom_ids )
        	 && ( distance_squared(
        		  $atom_site->{$atom_id},
        		  $atom_site->{$neighbour_id} ) <= $distance ** 2 ) ) {
        	    $around_atom_site{$neighbour_id} =
        		$atom_site->{$neighbour_id};
        	}
            }
        }
    }

    return \%around_atom_site;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Then, all atoms' distances are compared pairwisely in 1 + 26 surrounding
# boxes. If distance correspond to appropriate length, then connection is
# made by two atoms.
# Input:
#     $atom_site - atom data structure.
# Output:
#     connects atoms by adding "connection" key and values to atom site data
#     structure.
#

sub connect_atoms
{
    my ( $atom_site ) = @_;

    # Removes all previously described connections.
    for my $atom_id ( keys %{ $atom_site } ) {
	if( exists $atom_site->{$atom_id}{'connections'} ) {
	    delete $atom_site->{$atom_id}{'connections'}
	}
    }

    # For each cell, checks neighbouring cells. Creates box around atoms, makes
    # grid with edge length of max covalent radii of the parameter file.
    my ( $grid_box, undef ) = grid_box( $atom_site );
    my $neighbour_cells = identify_neighbour_cells( $grid_box );

    foreach my $cell ( keys %{ $grid_box } ) {
        foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
            foreach my $neighbour_id ( @{ $neighbour_cells->{$cell} } ) {
                if( ( is_connected( $atom_site->{"$atom_id"},
                                    $atom_site->{"$neighbour_id"} ) )
                 && ( ( ! exists $atom_site->{$atom_id}{'connections'} )
                   || ( ! any { $neighbour_id eq $_ }
                             @{ $atom_site->{$atom_id}{'connections'} } ) ) ){
                    push( @{ $atom_site->{$atom_id}{'connections'} },
                          "$neighbour_id" );
                }
            }
        }
    }

    return;
}

1;
