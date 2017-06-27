package ConnectAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( connect_atoms
                     create_box );

use List::Util qw( max min );
use Data::Dumper;

use lib qw( ./ );
use CifParser qw( select_atom_data );
use Combinatorics qw( permutation );
use LoadParams qw( covalent_radii );

my $parameter_file = "../../parameters/covalent_radii.csv";

# ------------------------------ Connect atoms ------------------------------- #

#
# Parameters
#

my %COVALENT_RADII = %{ covalent_radii( $parameter_file ) };

my $MAX_BOND_LENGTH =
    max( map { $COVALENT_RADII{$_}{"bond_length"} }
	 keys( %COVALENT_RADII ) ) * 2;

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
# Then, all atoms' distances are compared pairwisely in one box + 26 surrounding
# boxes. If distance is correspond to appropriate length, then connection is
# made by two atoms.
# Input  (1 arg): cif data structure form (look at CifParser.pm).
# Output (1 arg): atom data in cif data structure form that has additional
#                 data for each atom - hash of atom coordinates (x, y, z) as
#                 keys and atom coordinates that are connected to as values.
#

sub connect_atoms
{
    my ( $atom_site ) = @_;

    my $all_atom_coord =
    	select_atom_data( [ "id", "Cartn_x", "Cartn_y", "Cartn_z" ],
    			  $atom_site );

    # Creates smallest box that contain all atoms.
    my @boundary_box =
	create_box( @{ select_atom_data( [ "Cartn_x", "Cartn_y", "Cartn_z" ],
					 $atom_site ) } );

    # Assign atoms to cells in grid box.
    my %grid_box; # Box with cells with edge length of largest covalent radii
                  # in parameter file.
    my $cell_index_x;
    my $cell_index_y;
    my $cell_index_z;

    foreach my $atom_coord ( @{ $all_atom_coord } ) {
    	$cell_index_x =
    	    int( ( $atom_coord->[1] - $boundary_box[0] )
		 / $MAX_BOND_LENGTH ) + 1;
    	$cell_index_y =
    	    int( ( $atom_coord->[2] - $boundary_box[2] )
		 / $MAX_BOND_LENGTH ) + 1;
    	$cell_index_z =
    	    int( ( $atom_coord->[3] - $boundary_box[4] )
		 / $MAX_BOND_LENGTH ) + 1;

    	if( exists $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} ) {
    	    push( @{ $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} },
    		  $atom_coord );
    	} else {
    	          $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} =
    		  [ $atom_coord ];
    	}
    }

    # For each cell, checks neighbouring cells.
    my %connected_atoms = %{ $atom_site };
    my @cell_idx;

    foreach my $cell ( keys %grid_box ) {
    	@cell_idx = split( ",", $cell );
    	my @neighbour_cells; # The array will contain all atoms of the
                             # neighbouring 26 cells.

    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_idx[0] - 1..$cell_idx[0] + 1 ) ) {
    	for my $j ( ( $cell_idx[1] - 1..$cell_idx[1] + 1 ) ) {
    	for my $k ( ( $cell_idx[2] - 1..$cell_idx[2] + 1 ) ) {
	if( exists $grid_box{"$i,$j,$k"} ) {
	    push( @neighbour_cells, @{ $grid_box{"$i,$j,$k"} } ); } } } }

	# Checks, if there are connections between atoms.
	my $atom_type;
	my $neighbour_type;

	my $bond_length_comb;
	my $bond_length;
	my $length_error_comb;
	my $length_error;
	my $distance_btw_atoms;

        # Atoms that have been already checked for connections.
	my @checked_atoms;

    	foreach my $atom_coord ( @{ $grid_box{$cell} } ) {
    	    $connected_atoms{"data"}{$atom_coord->[0]}{"connections"} = [];
	    push( @checked_atoms, $atom_coord->[0] ); # Marks as visited atom.
	    $atom_type =
		$atom_site->{"data"}{$atom_coord->[0]}{"type_symbol"};
    	    foreach my $neighbour_atom ( @neighbour_cells ) {
    		# Checks distance between neighbouring atoms by formula:
    		# x^2+y^2+z^2 < (bond_length)^2.
		$neighbour_type =
		    $atom_site->{"data"}{$neighbour_atom->[0]}{"type_symbol"};
		$bond_length_comb =
		    permutation(
			2,
			[],
			[ $COVALENT_RADII{$atom_type}{"bond_length"},
			  $COVALENT_RADII{$neighbour_type}{"bond_length"} ],
			[] );
		$length_error_comb =
		    permutation(
			2,
			[],
			[ $COVALENT_RADII{$atom_type}{"length_error"},
			  $COVALENT_RADII{$neighbour_type}{"length_error"} ],
			[] );

    	    	$distance_btw_atoms =
    		    ( $neighbour_atom->[1] - $atom_coord->[1] ) ** 2
    	    	  + ( $neighbour_atom->[2] - $atom_coord->[2] ) ** 2
    	    	  + ( $neighbour_atom->[3] - $atom_coord->[3] ) ** 2;

    		for( my $i = 0; $i < scalar( @{ $bond_length_comb } ); $i++ ) {
    		    $bond_length =
			$bond_length_comb->[$i][0]
		      + $bond_length_comb->[$i][1];
    		    $length_error =
			$length_error_comb->[$i][0]
		      + $length_error_comb->[$i][1];
    		    if( ( $distance_btw_atoms >
    		          ( $bond_length - $length_error ) ** 2 )
    		     && ( $distance_btw_atoms <
    		          ( $bond_length + $length_error ) ** 2 ) ) {
    		        push( @{ $connected_atoms{"data"}
    		    		 {$atom_coord->[0]}
    		    		 {"connections"} },
    		    	      $neighbour_atom->[0] );
    		    	last;
    		    }
    		}
    	    }

	    my $cell_atom_type;

	    # Chekcs for connections inside current cell in grid box.
	    foreach my $cell_atom_coord ( @{ $grid_box{$cell} } ) {
		if( $atom_coord->[0] != $cell_atom_coord->[0]
		    && ! grep { $atom_coord->[0] } @checked_atoms ) {
		    $cell_atom_type =
			$atom_site->{"data"}
		                    {$cell_atom_coord->[0]}
		                    {"type_symbol"};

		    $bond_length_comb =
			permutation(
			    2,
			    [],
			    [ $COVALENT_RADII{$atom_type}{"bond_length"},
			      $COVALENT_RADII{$cell_atom_type}{"bond_length"} ],
			    [] );
		    $length_error_comb =
			permutation(
			    2,
			    [],
			    [ $COVALENT_RADII{$atom_type}{"length_error"},
			      $COVALENT_RADII{$cell_atom_type}{"length_error"} ],
			    [] );

		    $distance_btw_atoms =
			  ( $cell_atom_coord->[1] - $atom_coord->[1] ) ** 2
			+ ( $cell_atom_coord->[2] - $atom_coord->[2] ) ** 2
			+ ( $cell_atom_coord->[3] - $atom_coord->[3] ) ** 2;

		    for( my $i = 0; $i < scalar( @{ $bond_length_comb } ); $i++ ) {
			$bond_length =
			    $bond_length_comb->[$i][0]
			  + $bond_length_comb->[$i][1];
			$length_error =
			    $length_error_comb->[$i][0]
			  + $length_error_comb->[$i][1];
			if( ( $distance_btw_atoms >
			      ( $bond_length - $length_error ) ** 2 )
			    && ( $distance_btw_atoms <
				 ( $bond_length + $length_error ) ** 2 ) ) {
			    push( @{ $connected_atoms{"data"}
				     {$atom_coord->[0]}
				     {"connections"} },
				  $cell_atom_coord->[0] );
			    last;
			}
		    }
		}
	    }
    	}
    }

    return \%connected_atoms;
}

1;
