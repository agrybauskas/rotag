package ConnectAtoms;

use strict;
use warnings;

use List::Util qw(min max);

use lib qw( ./ );
use CifParser;

# ------------------------------ Connect atoms ------------------------------- #

#
# Shows what atom is connected to what atom using only information about atom
# coordinates.
#

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
# Input  (1 arg): atom coordinates in x, y, z form.
# Output (6 arg): coordinates of min and max x, y, z box boundaries in which
#                  all given atoms are contained.
#

sub create_box
{
    my @atom_coord = @_;

    my @atom_coord_x = map { $_->[0] } @atom_coord;
    my @atom_coord_y = map { $_->[1] } @atom_coord;
    my @atom_coord_z = map { $_->[2] } @atom_coord;

    # Directions are adapted to right-handed Cartesian coordinate system.
    # Looking for leftmost and rightmost coordinates of X-axis.
    my $most_left_x_coord     = min( @atom_coord_x );
    my $most_right_x_coord    = max( @atom_coord_x );

    # Looking for most backward and forward coordinates of Y-axis.
    my $most_backward_y_coord = min( @atom_coord_y );
    my $most_forward_y_coord  = max( @atom_coord_y );

    # Looking for downmost and upmost coordinates of Z-axis.
    my $most_down_z_coord     = min( @atom_coord_z );
    my $most_up_z_coord       = max( @atom_coord_z );

    # Coordinates of minimum bounding box that contains all given atoms.
    return $most_left_x_coord,     $most_right_x_coord,
           $most_backward_y_coord, $most_forward_y_coord,
           $most_down_z_coord,     $most_up_z_coord;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Then, all atoms' distances are compared pairwisely in one box. If distance
# is correspond to appropriate length, then connection is made by two atoms.
# Input  (2 arg): bond length in angstroms, bond length error and atom data
#                 in cif data structure form (look at CifParser.pm).
# Output (1 arg): atom data in cif data structure form that has additional
#                 data for each atom - hash of atom coordinates (x, y, z) as
#                 keys and atom coordinates that are connected to as values.
#

sub connect_atoms
{
    # TODO: bond length and error must be adjusted to variety of atoms.
    my ( $bond_length, $length_error, $atom_site ) = @_;

    my $all_atom_coord =
	CifParser::select_atom_data( [ "id", "Cartn_x", "Cartn_y", "Cartn_z" ],
				     $atom_site );

    # Creates smallest box that contain all atoms.
    my $only_atom_coord =
	CifParser::select_atom_data( [ "Cartn_x", "Cartn_y", "Cartn_z" ],
				     $atom_site );
    my @boundary_box = create_box( @$only_atom_coord );

    my %grid_box;

    my $cell_index_x;
    my $cell_index_y;
    my $cell_index_z;

    my %connected_atoms = %$atom_site;

    # Assign atoms to cells in grid box.
    foreach my $atom_coord ( @$all_atom_coord ) {
    	$cell_index_x =
    	    int( ( $atom_coord->[1] - $boundary_box[0] ) / $bond_length ) + 1;
    	$cell_index_y =
    	    int( ( $atom_coord->[2] - $boundary_box[2] ) / $bond_length ) + 1;
    	$cell_index_z =
    	    int( ( $atom_coord->[3] - $boundary_box[4] ) / $bond_length ) + 1;

    	if( exists $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} ) {
    	    push( @{ $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} },
    		  $atom_coord );
    	} else {
    	          $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} =
    		  [ $atom_coord ];
    	}
    }

    my @cell_idx;

    # For each cell, checks neighbouring cells.
    foreach my $cell ( keys %grid_box ) {
    	@cell_idx = split( ",", $cell );

    	my @neighbour_cells; # The array will contain all atoms of the
                             # neighbouring 26 cells.

    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_idx[0] - 1..$cell_idx[0] + 1 ) ) {
    	    for my $j ( ( $cell_idx[1] - 1..$cell_idx[1] + 1 ) ) {
    		for my $k ( ( $cell_idx[2] - 1..$cell_idx[2] + 1 ) ) {
    		    if( exists $grid_box{"$i,$j,$k"} ) {
    			push( @neighbour_cells, @{ $grid_box{"$i,$j,$k"} } );
    		    }
    		}
    	    }
    	}

    	# TODO: add atoms that are in the center cell.
    	foreach my $cell_atom_coord ( @{ $grid_box{$cell} } ) {
	    $connected_atoms{"data"}{$cell_atom_coord->[0]}{"connections"} = [];
    	    foreach my $neighbour_atom ( @neighbour_cells ) {
    		# Checks distance between neighbouring atoms by formula:
    		# x^2+y^2+z^2 < (bond_length)^2
    		my $distance_btw_atoms;

    	    	$distance_btw_atoms =
    		    ( $neighbour_atom->[1] - $cell_atom_coord->[1] ) ** 2
    	    	  + ( $neighbour_atom->[2] - $cell_atom_coord->[2] ) ** 2
    	    	  + ( $neighbour_atom->[3] - $cell_atom_coord->[3] ) ** 2;
    		if( ( $distance_btw_atoms >
    		      ( $bond_length - $length_error ) ** 2 )
    		 && ( $distance_btw_atoms <
    		      ( $bond_length + $length_error ) ** 2 ) ) {
		    push( @{ $connected_atoms{"data"}
			                     {$cell_atom_coord->[0]}
			                     {"connections"} },
			  $neighbour_atom->[0] );
    		}
    	    }
    	}
    }

    return \%connected_atoms;
}

1;
