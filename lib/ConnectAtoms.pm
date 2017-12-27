package ConnectAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( bond_type
                     connect_atoms
                     create_box
                     grid_box
                     hybridization
                     is_connected
                     is_second_neighbour );

use List::Util qw( any
                   max
                   min );

use AtomProperties qw( %ATOMS );
use Combinatorics qw( permutation );
use PDBxParser qw( select_atom_data );
use MoleculeProperties qw( %BOND_TYPES );

# ------------------------------ Connect atoms ------------------------------- #

#
# Shows what atom is connected to what atom using only information about atom
# coordinates.
#

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
# Input:
#     @atom_coord - list of atom coordinates in x, y, z array form.
# Output:
#     coordinates of min and max x, y, z box boundaries in which all given atoms
#     are contained.
#

sub create_box
{
    my ( $atom_coord ) = @_;

    my @atom_coord_x = map { $_->[0] } @{ $atom_coord };
    my @atom_coord_y = map { $_->[1] } @{ $atom_coord };
    my @atom_coord_z = map { $_->[2] } @{ $atom_coord };

    # Directions are adapted to right-handed Cartesian coordinate system.
    # Looking for leftmost and rightmost coordinates of X-axis.
    my $min_coord_x  = min( @atom_coord_x );
    my $max_coord_x  = max( @atom_coord_x );

    # Looking for most backward and forward coordinates of Y-axis.
    my $min_coord_y  = min( @atom_coord_y );
    my $max_coord_y  = max( @atom_coord_y );

    # Looking for downmost and upmost coordinates of Z-axis.
    my $min_coord_z  = min( @atom_coord_z );
    my $max_coord_z  = max( @atom_coord_z );

    # Coordinates of minimum bounding box that contains all given atoms.
    return [ $min_coord_x, $max_coord_x,
	     $min_coord_y, $max_coord_y,
	     $min_coord_z, $max_coord_z ];
}

#
# Divides atoms into grid box of given edge length.
# Input:
#     $atom_site - special data structure.
#     $edge_length - edge length of the cell inside grid box.
# Output:
#     %grid_box - hash where key is string representing cell id and value -
#     atom id.
#

sub grid_box
{
    my ( $atom_site, $edge_length ) = @_;

    # Default value for edge length is two times greater than the largest
    # covalent radius.
    $edge_length //=
	max( map { @{ $ATOMS{$_}{"covalent_radius"}{"length"} } }
	     keys %ATOMS ) * 2;

    # Determines boundary box around all atoms.
    my $atom_data =
	select_atom_data( $atom_site,
			  [ "id", "Cartn_x", "Cartn_y", "Cartn_z" ] );
    my @atom_coordinates = map { [ $_->[1], $_->[2], $_->[3] ] } @{ $atom_data };
    my $boundary_box = create_box( \@atom_coordinates );

    # Creates box with cells with edge length of given variable in angstroms.
    my %grid_box;
    my $cell_index_x;
    my $cell_index_y;
    my $cell_index_z;

    # Iterates through atoms and determines in which cell these atoms are.
    foreach my $atom_coord ( @{ $atom_data } ) {
    	$cell_index_x =
    	    int( ( $atom_coord->[1] - $boundary_box->[0] )
    		 / $edge_length ) + 1;
    	$cell_index_y =
    	    int( ( $atom_coord->[2] - $boundary_box->[2] )
    		 / $edge_length ) + 1;
    	$cell_index_z =
    	    int( ( $atom_coord->[3] - $boundary_box->[4] )
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

#
# Checks if two atoms are connected.
# Input:
#     $target_atom - atom data structure (see PDBxParser) of first atom.
#     $neighbour_atom - atom data structure of second atom.
# Output:
#     $is_connected - boolean of two values: 0 (as false) and 1 (as true).
#

sub is_connected
{
    my ( $target_atom, $neighbour_atom ) = @_;

    # Generates all possible combinations of covalent distances.
    my $bond_length_comb =
    	permutation( 2,
    		     [],
    		     [ $ATOMS{$target_atom->{"type_symbol"}}
		             {"covalent_radius"}
		             {"length"},
		       $ATOMS{$neighbour_atom->{"type_symbol"}}
		             {"covalent_radius"}
		             {"length"} ],
    		     [] );

    my $length_error_comb =
    	permutation( 2,
    		     [],
    		     [ $ATOMS{$target_atom->{"type_symbol"}}
		             {"covalent_radius"}
		             {"error"},
		       $ATOMS{$neighbour_atom->{"type_symbol"}}
		             {"covalent_radius"}
		             {"error"} ],
    		     [] );

    # Precalculates distance between atom pairs.
    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

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
	@{ $atom_site->{"$target_atom_id"}{"connections"} } ) {
    foreach my $j (
	@{ $atom_site->{$i}{"connections"} } ) {
	if( "$sec_neighbour_id" eq "$j" ) {
	    $is_sec_neighbour = 1;
	    last;
	}
    } last if $is_sec_neighbour == 1; }

    return $is_sec_neighbour;
}

sub bond_type
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $target_atom_type = $target_atom->{"type_symbol"};
    my $neighbour_atom_type = $neighbour_atom->{"type_symbol"};

    # Precalculates squared distance between atom pairs. Delocalized bonds are
    # described by double or triple bond.
    # TODO: investigate, if this delocalized bond simplification can be made.
    my $squared_distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    for my $bond_type ( keys %BOND_TYPES ) {
	if( exists $BOND_TYPES{$bond_type}
	                      {$target_atom_type}
	                      {$neighbour_atom_type}
	 || exists $BOND_TYPES{$bond_type}
	                      {$neighbour_atom_type}
	                      {$target_atom_type} ) {
	    my $bond_length_min = $BOND_TYPES{$bond_type}
	                                     {$target_atom_type}
                                	     {$neighbour_atom_type}
    	                                     {"min_length"}
    	                       || $BOND_TYPES{$bond_type}
                                 	     {$neighbour_atom_type}
                                	     {$target_atom_type}
    	                                     {"min_length"};
	    my $bond_length_max = $BOND_TYPES{$bond_type}
                                 	     {$target_atom_type}
                                	     {$neighbour_atom_type}
    	                                     {"max_length"}
    	                       || $BOND_TYPES{$bond_type}
                                 	     {$neighbour_atom_type}
                                	     {$target_atom_type}
    	                                     {"max_length"};

	    if( ( $squared_distance >  $bond_length_min**2 )
	     && ( $squared_distance <= $bond_length_max**2 ) ) {
	    	return $bond_type;
	    	last;
	    }
	}
    }
}

sub hybridization
{
    # Use connect_atoms before using hybridization function.
    my ( $atom_site ) = @_;

    for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
	# Determines every type of connection.
	my @bond_types;
	for my $connection_id ( @{ $atom_site->{$atom_id}{"connections"} } ) {
	    push( @bond_types,
		  bond_type( $atom_site->{$atom_id},
			     $atom_site->{$connection_id} ) );
	}

	# Depending on connections, assigns hybridization type.
	# TODO: check more possibilities of different bonds and their
	# combinations.
	if( any { $_ eq "double" } @bond_types ) {
	    $atom_site->{$atom_id}{"hybridization"} = "sp2";
	} elsif( any { $_ eq "triple" } @bond_types ) {
	    $atom_site->{$atom_id}{"hybridization"} = "sp";
	} else {
	    $atom_site->{$atom_id}{"hybridization"} = "sp3";
	}
    }

    return $atom_site;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Then, all atoms' distances are compared pairwisely in one box + 26 surrounding
# boxes. If distance is correspond to appropriate length, then connection is
# made by two atoms.
# Input:
#     $atom_site - atom data structure.
# Output:
#     %connected_atoms - atom data in pbdx data structure form that has
#     additional data for each atom - hash of atom coordinates (x, y, z) as keys
#     and atom coordinates that are connected to as values.
#

sub connect_atoms
{
    my ( $atom_site ) = @_;

    # For each cell, checks neighbouring cells.
    my %connected_atoms = %{ $atom_site };
    my @cell_indexes;

    # Creates box around atoms, makes grid with edge length of max covalent radii
    # of the parameter file.
    my $grid_box = grid_box( $atom_site );

    # Checks for neighbouring cells for each cell.
    foreach my $cell ( keys %{ $grid_box } ) {
    	@cell_indexes = split( ",", $cell );
	my @neighbour_cells; # The array will contain all atoms of the
	                     # neighbouring 26 cells.
    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_indexes[0] - 1..$cell_indexes[0] + 1 ) ) {
    	for my $j ( ( $cell_indexes[1] - 1..$cell_indexes[1] + 1 ) ) {
    	for my $k ( ( $cell_indexes[2] - 1..$cell_indexes[2] + 1 ) ) {
    	if( exists $grid_box->{"$i,$j,$k"} ) {
    	    push( @neighbour_cells, @{ $grid_box->{"$i,$j,$k"} } ); } } } }

    	# Checks, if there are connections between atoms.
    	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
    	    foreach my $neighbour_id ( @neighbour_cells ) {
    		if( is_connected( $atom_site->{"$atom_id"},
    				  $atom_site->{"$neighbour_id"} ) ) {
    		    push( @{ $connected_atoms{$atom_id}{"connections"} },
    		          "$neighbour_id" );
    		}
    	    }
    	}
    }

    return \%connected_atoms;
}

1;
