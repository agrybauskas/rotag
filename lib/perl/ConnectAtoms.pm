package ConnectAtoms;

use strict;
use warnings;

use List::Util qw(min max);
use Data::Dumper;

use Math::Round;

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
# Input   (1 arg): atom coordinates in x, y, z form.
# Output (6 args): coordinates of min and max x, y, z box boundaries in which
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
# Input (2 args): bond length in angstroms, coordinates of atoms (x, y, z).
#

sub connect_atoms
{
    my $bond_length = shift;
    my @all_atom_coord  = @_;

    # Creates smallest box that contain all atoms.
    my @boundary_box = create_box( @all_atom_coord );

    # Divides box into cells, edges of which has length of bond.
    my $cell_number_x =
    	int( ( $boundary_box[1] - $boundary_box[0] ) / $bond_length ) + 1;
    my $cell_number_y =
    	int( ( $boundary_box[3] - $boundary_box[2] ) / $bond_length ) + 1;
    my $cell_number_z =
    	int( ( $boundary_box[5] - $boundary_box[4] ) / $bond_length ) + 1;

    my @grid_box;

    foreach my $idx_x ( ( 1..$cell_number_x ) ) {
    	foreach my $idx_y ( ( 1..$cell_number_y ) ) {
    	    foreach my $idx_z ( ( 1..$cell_number_z ) ) {
    		push( @grid_box, [ $idx_x, $idx_y, $idx_z ], [] );
    	    }
    	}
    }

    # Assigns atoms to cells in grid box.
    # foreach my $atom_coord ( @all_atom_coord ) {
    # 	for( my $cell = 0; $cell < $#grid_box; $cell += 2 ) {
    # 	    if( $atom_coord->[0] <
    # 		$boundary_box[0] + ( $grid_box[$cell]->[0] - 1 ) * $bond_length ) {
    # 		next;
    # 	    }

    # 	    if( $atom_coord->[0] >=
    # 		$boundary_box[0] +   $grid_box[$cell]->[0] * $bond_length ) {
    # 		next;
    # 	    }

    # 	    if( $atom_coord->[1] <
    # 		$boundary_box[2] + ( $grid_box[$cell]->[1] - 1 ) * $bond_length ) {
    # 		next;
    # 	    }

    # 	    if( $atom_coord->[1] >=
    # 		$boundary_box[2] +   $grid_box[$cell]->[1] * $bond_length ) {
    # 		next;
    # 	    }

    # 	    if(	$atom_coord->[2] <
    # 		$boundary_box[4] + ( $grid_box[$cell]->[2] - 1 ) * $bond_length ) {
    # 		next;
    # 	    }

    # 	    if( $atom_coord->[2] >=
    # 		$boundary_box[4] +   $grid_box[$cell]->[2] * $bond_length ) {
    # 		push( @{ $grid_box[$cell+1] }, [ @$atom_coord ] );
    # 	    }
    # 	}
    # }

    # Removes empty cells.
    # my @grid_box_full; # Without empty cells.
    # my %grid_box_full; # Hash version of the same array. Will be esier to find
    #                    # neighbouring cells.

    # for( my $cell = 0; $cell < $#grid_box; $cell += 2 ) {
    #     if( $#{ $grid_box[$cell+1] } != -1 ) {
    #         push( @grid_box_full, 
    # 		  join( ",", @{ $grid_box[$cell] } ), ( $grid_box[$cell+1] ) );
    #     }
    # }

    # %grid_box_full = @grid_box_full;

    # my @cell_idx;

    # For each cell, checks neighbouring cells.
    # foreach my $cell ( keys %grid_box_full ) {
    # 	@cell_idx = split( ",", $cell );

    # 	my @neighbour_cells; # The array will contain all atoms of the neighbouring
    # 	                     # 26 cells.

	# TODO: maybe introduce regular expression instead of looking for every
	# possible cell.
	# $i represents x, $j - y, $k - z coordinates.
	# for my $i ( ( $cell_idx[0] - 1..$cell_idx[0] + 1 ) ) {
	#     for my $j ( ( $cell_idx[1] - 1..$cell_idx[1] + 1 ) ) {
	# 	for my $k ( ( $cell_idx[2] - 1..$cell_idx[2] + 1 ) ) {      
	# 	    if ( exists $grid_box_full{"$i,$j,$k"} ) {
	# 		push( @neighbour_cells, @{ $grid_box_full{"$i,$j,$k"} } );
	# 	    }
	# 	}
	#     }
	# }

	# Checks distance between neighbouring atoms by formula:
	# x^2+y^2+z^2 < (bond_length)^2
	# my $distance_btw_atoms;

	# TODO: add atoms that are in the center cell.
    	# foreach my $cell_atom_coord ( @{ $grid_box_full{$cell} } ) {
    	#     foreach my $neighbour_atom ( @neighbour_cells ) {
	# 	$distance_btw_atoms =  
	# 	      ( $neighbour_atom->[0] - $cell_atom_coord->[0] ) ** 2
	# 	    + ( $neighbour_atom->[1] - $cell_atom_coord->[1] ) ** 2
	# 	    + ( $neighbour_atom->[2] - $cell_atom_coord->[2] ) ** 2;

	# 	# TODO: deal with approximate distances in mmCif.
	# 	if( ( $distance_btw_atoms >
	# 	      ( $bond_length - 0.1 * $bond_length ) ** 2 )
	# 	 && ( $distance_btw_atoms <
	# 	      ( $bond_length + 0.1 * $bond_length ) ** 2 ) ) {
	# 	    print sqrt( $distance_btw_atoms ) . "\n";
	# 	}
    	
	#     }
	# }
    # }
}

1;
