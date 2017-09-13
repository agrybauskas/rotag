package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( any max );

use lib qw( ./ );
use PDBxParser qw( filter_atoms select_atom_data );
use ConnectAtoms qw( connect_atoms grid_box is_connected is_second_neighbour );
use LoadParams qw( rotatable_bonds vdw_radii );

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

my $MAX_VDW_RADIUS =
    max( map { vdw_radii()->{$_} }
	 keys( %{ vdw_radii() } ) ) * 2;

#
# Checks, if two atoms are colliding.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#          $is_colliding - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_colliding
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $vdw_length =
	vdw_radii()->{$target_atom->{"type_symbol"}}
      + vdw_radii()->{$neighbour_atom->{"type_symbol"}};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    # Checks, if distance between atom pairs is in one of the combinations.
    my $is_colliding;

    if( $distance < $vdw_length ** 2 ) {
	$is_colliding = 1;
    } else {
	$is_colliding = 0;
    }

    return $is_colliding;
}

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
# Input:
#     $atom_site - atom data structure.
#     $atom_specifier - compound data structure for specifying desirable atoms
#     (see PDBxParser.pm).
# Output:
#     %atom_clashes - modified $atom_site with added information about atom
#     clashes of specified atoms.
#

sub radius_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Clashes of all atoms analyzed, if no specific atoms are selected.
    $atom_specifier = { "group_pdb" => [ "ATOM" ] } unless $atom_specifier;
    my @atom_ids = # Atom ids selected by $atom_specifier.
    	map { $_->[0] }
        @{ select_atom_data( filter_atoms( $atom_site, $atom_specifier ),
			     [ "id" ] ) };

    # Checks for connecting atoms that will be excluded from clash list.
    connect_atoms( $atom_site );

    # Creates box around atoms, makes grid with edge length of max covalent
    # radii.
    my $grid_box = grid_box( $atom_site, $MAX_VDW_RADIUS );

    # Checks for neighbouring cells for each cell.
    foreach my $cell ( keys %{ $grid_box } ) {
    	my @cell_indexes = split( ",", $cell );
    	my @neighbour_cells; # The array will contain all atoms of the
                             # neighbouring 26 cells.

    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_indexes[0] - 1..$cell_indexes[0] + 1 ) ) {
    	for my $j ( ( $cell_indexes[1] - 1..$cell_indexes[1] + 1 ) ) {
    	for my $k ( ( $cell_indexes[2] - 1..$cell_indexes[2] + 1 ) ) {
    	if( exists $grid_box->{"$i,$j,$k"} ) {
    	    push( @neighbour_cells, @{ $grid_box->{"$i,$j,$k"} } ); } } } }

    	# Checks, if there are clashes between atoms.
    	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
    	    if( any { $atom_id eq $_ } @atom_ids ) {
    		foreach my $neighbour_id ( @neighbour_cells ) {
    		    if( not any { $neighbour_id eq $_ } @atom_ids ) {
    		    	if( is_colliding( $atom_site->{"$atom_id"},
    		    			  $atom_site->{"$neighbour_id"} )
    		    	    && $atom_id ne $neighbour_id
    		    	    && ( not is_connected(
    				     $atom_site->{"$atom_id"},
    				     $atom_site->{"$neighbour_id"} ) )
    		    	    && ( not is_second_neighbour(
    				     $atom_site, $atom_id, $neighbour_id ) ) ) {
    		    	    push( @{ $atom_site->{$atom_id}{"clashes"} },
    		    		  $neighbour_id );
    		    	}
    		    }
    		}
    	    }
    	}
    }

    # Removes atoms with any clashes.
    foreach my $atom_id ( keys %{ $atom_site } ) {
    	delete $atom_site->{$atom_id}
    	if exists $atom_site->{$atom_id}{"clashes"};
    }

    return $atom_site;
}

1;
