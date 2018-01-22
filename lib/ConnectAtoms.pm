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
                     is_neighbour
                     is_second_neighbour
                     rotatable_bonds
                     sort_by_priority );

use List::MoreUtils qw( uniq );
use List::Util qw( any
                   max
                   min );

use AtomProperties qw( %ATOMS );
use Combinatorics qw( permutation );
use PDBxParser qw( filter );
use LinearAlgebra qw( flatten );
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
	filter( { "atom_site" => $atom_site,
		  "data" => [ "id", "Cartn_x", "Cartn_y", "Cartn_z" ] } );
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

sub is_neighbour
{
    my ( $atom_site, $target_atom_id, $neighbour_id ) = @_;

    my $is_neighbour = 0;
    foreach my $i ( @{ $atom_site->{"$target_atom_id"}{"connections"} } ) {
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
}

sub sort_by_priority
{
    my ( $atom_names ) = @_;

    # First priority is decided by atom type: S > P > O > N > C > H.
    # Second priority - by greek letter: A > B > G > D > E > Z > H.
    # TODO: look for more greek letters that are in PDBx.
    # Third priority - by numeration: 1 > 2 > 3 and etc.
    # This priority is achieved by assinging first, second and third priorities
    # to numbers. Then iteratively is sorted by priorities.
    my %atom_type_priority =
	( "H" => 1, "C" => 2, "N" => 3, "O" => 4, "P" => 5, "S" => 6 );
    my %greek_letter_priority =
	( "H" => 1, "Z" => 2, "E" => 3, "D" => 4, "G" => 5, "B" => 6,
	  "A" => 7,  "" => 8 );

    # Decomposes each atom name by its components.
    my %atom_names;
    for my $atom_name ( @{ $atom_names } ) {
	my ( $atom_type ) = $atom_name =~ /(^[a-zA-z])[a-zA-z]?\d?/;
	my ( $greek_letter ) = $atom_name =~ /^[a-zA-z]([a-zA-z]?)\d?/;
	my ( $number ) = $atom_name =~ /^[a-zA-z][a-zA-z]?(\d?)/;
	$atom_names{$atom_name}{"type"} =
	    $atom_type_priority{$atom_type};
	$atom_names{$atom_name}{"greek_letter"} =
	    $greek_letter_priority{$greek_letter};
	$atom_names{$atom_name}{"number"} = $number;
    }

    # Sorts by rules of described in %atom_names.
    my @sorted_names =
    	sort {
	    $atom_names{$b}{"type"} <=> $atom_names{$a}{"type"}
	 || $atom_names{$b}{"greek_letter"} <=> $atom_names{$a}{"greek_letter"}
         || $atom_names{$a}{"number"} cmp $atom_names{$b}{"number"}
        } @{ $atom_names };

    return \@sorted_names;
}

sub rotatable_bonds
{
    my ( $atom_site, $start_atom_id, $next_atom_id ) = @_;

    # By default, CA is starting atom and CB next.
    $start_atom_id //= filter( { "atom_site" => $atom_site,
				 "include" => { "label_atom_id" => [ "CA" ] },
				 "data" => [ "id" ],
				 "is_list" => 1 } )->[0];
    $next_atom_id //=  filter( { "atom_site" => $atom_site,
				 "include" => { "label_atom_id" => [ "CB" ] },
				 "data" => [ "id" ],
				 "is_list" => 1 } )->[0];

    my %atom_site = %{ $atom_site }; # Copy of the variable.
    my @atom_ids = keys %atom_site;
    my @visited_atom_ids = ( $start_atom_id );
    my @next_atom_ids = ( $next_atom_id );
    my %parent_atom_ids;

    my %rotatable_bonds;

    # Marks parent atom for next atom id.
    $parent_atom_ids{$next_atom_id} = $start_atom_id;

    # Connects and determines hybridization for each atom.
    connect_atoms( \%atom_site );
    hybridization( \%atom_site );

    # Exists if there are no atoms that is not already visited.
    while( scalar( @next_atom_ids ) != 0 ) {
    	# Iterates through every neighbouring atom if it was not visited
    	# before.
    	my @neighbour_atom_ids;
    	for my $atom_id ( @next_atom_ids ) {
    	    my $parent_atom_id = $parent_atom_ids{$atom_id};

    	    if( $atom_site{$parent_atom_id}{"hybridization"} eq "sp3"
    		|| $atom_site{$atom_id}{"hybridization"} eq "sp3" ) {
    		# If last visited atom was sp3, then rotatable bonds from
    		# previous atom are copied and the new one is appended.
    		push( @{ $rotatable_bonds{$atom_id} },
    		      [ $parent_atom_id, $atom_id ] );
    		unshift( @{ $rotatable_bonds{$atom_id} },
    			 @{ $rotatable_bonds{$parent_atom_id} } )
    		    if exists $rotatable_bonds{$parent_atom_id};
    	    } else {
    		# If last visited atom is sp2 or sp, inherits its rotatable
    		# bonds, because double or triple bonds do not rotate.
    		unshift( @{ $rotatable_bonds{$atom_id} },
    			 @{ $rotatable_bonds{$parent_atom_id} } )
    		    if exists $rotatable_bonds{$parent_atom_id};
    	    }

    	    # Marks visited atoms.
    	    push( @visited_atom_ids, $atom_id );

    	    # Marks neighbouring atoms.
    	    push( @neighbour_atom_ids,
    		  @{ $atom_site{$atom_id}{"connections"} } );

    	    # Marks parent atoms for each neighbouring atom.
    	    for my $neighbour_atom_id ( @neighbour_atom_ids ) {
    		$parent_atom_ids{$neighbour_atom_id} = $atom_id
    		    if ( ! grep { $neighbour_atom_id eq $_ }
    			   @visited_atom_ids )
    		    # HACK: this exception might produce unexpected results.
    		    && ( ! exists $parent_atom_ids{$neighbour_atom_id} );
    	    }
    	}

    	# Determines next atoms that should be visited.
    	@next_atom_ids = (); # Resets value for the new ones to be appended.
    	for my $neighbour_atom_id ( uniq @neighbour_atom_ids ) {
    	    if( ( ! grep { $neighbour_atom_id eq $_ } @visited_atom_ids )
    		&& ( any { $neighbour_atom_id eq $_ } @atom_ids ) ) {
    		push( @next_atom_ids, $neighbour_atom_id );
    	    }
    	}
    }

    # Removes bonds, if they have the id of the target atom. Also, remove ids,
    # which have no rotatable bonds after previous filtering.
    for my $atom_id ( keys %rotatable_bonds ) {
    	my $last_bond_idx = $#{ $rotatable_bonds{$atom_id} };
    	if( ( $atom_id == $rotatable_bonds{$atom_id}[$last_bond_idx][0]
    	   || $atom_id == $rotatable_bonds{$atom_id}[$last_bond_idx][1] ) ) {
    	    pop( @{ $rotatable_bonds{$atom_id} } );
    	}
    	if( ! @{ $rotatable_bonds{$atom_id} } ) {
    	    delete $rotatable_bonds{$atom_id};
    	}
    }

    return \%rotatable_bonds;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Then, all atoms' distances are compared pairwisely in one box + 26 surrounding
# boxes. If distance is correspond to appropriate length, then connection is
# made by two atoms.
# Input:
#     $atom_site - atom data structure.
#

sub connect_atoms
{
    my ( $atom_site ) = @_;

    # Removes all previously described connections.
    for my $atom_id ( keys %{ $atom_site } ) {
	delete $atom_site->{$atom_id}{"connections"}
	    if exists $atom_site->{$atom_id}{"connections"};
    }

    # For each cell, checks neighbouring cells. Creates box around atoms, makes
    # grid with edge length of max covalent radii of the parameter file.
    my @cell_indexes;
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
    		if( ( is_connected( $atom_site->{"$atom_id"},
				    $atom_site->{"$neighbour_id"} ) )
		 && ( ( ! exists $atom_site->{$atom_id}{"connections"} )
		   || ( ! any { $neighbour_id eq $_ }
			     @{ $atom_site->{$atom_id}{"connections"} } ) )){
		    push( @{ $atom_site->{$atom_id}{"connections"} },
			  "$neighbour_id" );
    		}
    	    }
    	}
    }
}

1;
