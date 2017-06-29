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

my $covalent_file = "../../parameters/covalent_radii.csv";

# ------------------------------ Connect atoms ------------------------------- #

#
# Parameters
#

my %COVALENT_RADII = %{ covalent_radii( $covalent_file ) };

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

sub grid_box
{
    my ( $atom_site, $edge_length ) = @_;

    # Determines boundary box around all atoms.
    my $all_atom_coord =
	select_atom_data( [ "id", "Cartn_x", "Cartn_y", "Cartn_z" ],
			  $atom_site );
    my @boundary_box = create_box( $all_atom_coord );

    # Creates box with cells with edge length of given variable in angstroms.
    my %grid_box;
    my $cell_index_x;
    my $cell_index_y;
    my $cell_index_z;

    # Iterates through atoms and determines in which cell these atoms are.
    foreach my $atom_coord ( @{ $all_atom_coord } ) {
	$cell_index_x =
	    int( ( $atom_coord->[1] - $boundary_box[0] )
		 / $edge_length ) + 1;
	$cell_index_y =
	    int( ( $atom_coord->[2] - $boundary_box[2] )
		 / $edge_length ) + 1;
	$cell_index_z =
	    int( ( $atom_coord->[3] - $boundary_box[4] )
		 / $edge_length ) + 1;

	# Checks if hash keys already  exist.
	if( exists $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} ) {
	    push( @{ $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} },
		  $atom_coord->[0] );
	} else {
	    $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} =
		[ $atom_coord->[0] ];
	}
    }

    return \%grid_box;
}

sub check_distance
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $bond_length_comb =
    	permutation( 2,
    		     [],
    		     [ $COVALENT_RADII{$target_atom->{"type_symbol"}}
    		       {"bond_length"},
    		       $COVALENT_RADII{$neighbour_atom->{"type_symbol"}}
    		                      {"bond_length"}],
    		     [] );

    my $length_error_comb =
    	permutation( 2,
    		     [],
    		     [ $COVALENT_RADII{$target_atom->{"type_symbol"}}
    		                      {"length_error"},
    		       $COVALENT_RADII{$neighbour_atom->{"type_symbol"}}
    	                              {"length_error"}],
    		     [] );

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    my $bond_length;
    my $length_error;

    for( my $i = 0; $i < scalar( @{ $bond_length_comb } ); $i++ ) {
    	$bond_length =
    	    $bond_length_comb->[$i][0]
    	  + $bond_length_comb->[$i][1];
    	$length_error =
    	    $length_error_comb->[$i][0]
           + $length_error_comb->[$i][1];
	if( ( $distance >= ( $bond_length - $length_error ) ** 2 )
	 && ( $distance <= ( $bond_length + $length_error ) ** 2 ) ) {
	    return "connected";
	}
    }
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

    # For each cell, checks neighbouring cells.
    my %connected_atoms = %{ $atom_site };
    my @cell_idx;

    # Creates box around atoms, makes grid with edge length of max covalent radii
    # of the parameter file.
    my $grid_box = grid_box( $atom_site, $MAX_BOND_LENGTH );

    # Checks for neighbouring cells for each cell.
    foreach my $cell ( keys %{ $grid_box } ) {
    	@cell_idx = split( ",", $cell );
    	my @neighbour_cells; # The array will contain all atoms of the
                             # neighbouring 26 cells.

    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_idx[0] - 1..$cell_idx[0] + 1 ) ) {
    	for my $j ( ( $cell_idx[1] - 1..$cell_idx[1] + 1 ) ) {
    	for my $k ( ( $cell_idx[2] - 1..$cell_idx[2] + 1 ) ) {
	if( exists $grid_box->{"$i,$j,$k"} ) {
	    push( @neighbour_cells, @{ $grid_box->{"$i,$j,$k"} } ); } } } }

        # Atoms that have been already checked for connections.
	my @checked_atoms;

	# Checks, if there are connections between atoms.
    	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
	    push( @checked_atoms, $atom_id ); # Marks as visited atom.
    	    foreach my $neighbour_id ( @neighbour_cells ) {
		if( check_distance(
			$atom_site->{"data"}{"$atom_id"},
			$atom_site->{"data"}{"$neighbour_id"} )
		    eq "connected" ) {
		    push( @{ $connected_atoms{"data"}
			                     {$atom_id}
			                     {"connections"} },
		          $neighbour_id );
		}
    	    }
    	}
    }

    return \%connected_atoms;
}

1;
