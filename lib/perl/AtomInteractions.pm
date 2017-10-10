package AtomInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( potential );

use List::Util qw( any max );

use lib qw( ./ );
use AtomProperties qw( %ATOMS );
use ConnectAtoms qw( connect_atoms
                     grid_box
                     is_connected
                     is_second_neighbour );
use PDBxParser qw( filter_atoms
                   select_atom_data );

# ---------------------------- General potential ------------------------------ #

#
# General potential function that inherits other potential functions. Its purpose
# is to calculate attractive/repulsive forces.
# Input:
#     $atom_site - atom data structure.
#     $potential - potential function that is used in calculation of forces.
#     $cutoff - value of potential where placement of pseudo atom is not
#     considered.
#     $target_atoms, $visible_atoms - compound data structure for specifying
#     desirable atoms (see PDBxParser.pm).
# Output:
#     %atom_site - modified $atom_site with added information about atom
#     energy value.
#

sub potential
{
    my ( $atom_site, $potential, $cutoff, $target_atoms, $visible_atoms ) = @_;

    # Interactions of all atoms analyzed, if no specific atoms are selected.
    $target_atoms =  { "group_pdb" => [ "ATOM" ] } unless $target_atoms;
    $visible_atoms = { "group_pdb" => [ "ATOM" ] } unless $visible_atoms;
    my @target_atom_ids = # Atom ids selected by $atom_specifier.
    	map { $_->[0] }
        @{ select_atom_data(
	   filter_atoms( $atom_site, $target_atoms ),
	   [ "id" ] ) };
    my @visible_atom_ids = # Atom ids selected by $atom_specifier.
    	map { $_->[0] }
        @{ select_atom_data(
	   filter_atoms( $atom_site, $visible_atoms ),
	   [ "id" ] ) };

    # Selection of potential function.
    my $potential_function;
    $potential_function = \&hard_sphere   if $potential eq "hard_sphere";
    $potential_function = \&soft_sphere   if $potential eq "soft_sphere";
    $potential_function = \&exponential   if $potential eq "exponential";
    $potential_function = \&leonard_jones if $potential eq "leonard_jones";

    # Checks for connecting atoms that will be excluded from clash list.
    connect_atoms( $atom_site );

    # Makes grid with edge length of max covalent radii.
    my $grid_box = grid_box( $atom_site );

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

    	# Calculates pair interactions of two atoms. First and second neighbours
    	# are not included in the analysis.
    	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
    	    $atom_site->{$atom_id}{"potential_energy"} = 0;
    	    if( any { $atom_id eq $_ } @target_atom_ids ) {
    	foreach my $neighbour_id ( @neighbour_cells ) {
    	    if( any { $neighbour_id eq $_ } @visible_atom_ids ) {
    	    if( $atom_id ne $neighbour_id
    		&& ( not is_connected( $atom_site->{"$atom_id"},
    				       $atom_site->{"$neighbour_id"} ) )
    		&& ( not is_second_neighbour( $atom_site,
    					      $atom_id,
    					      $neighbour_id ) ) ) {
    		$atom_site->{$atom_id}{"potential_energy"} +=
    		    $potential_function->( $atom_site->{"$atom_id"},
    					   $atom_site->{"$neighbour_id"} );
    		last if $atom_site->{$atom_id}{"potential_energy"} > $cutoff;
    	} } } } } }

    # Removes atoms with any clashes.
    foreach my $atom_id ( keys %{ $atom_site } ) {
    	delete $atom_site->{$atom_id}
    	if exists $atom_site->{$atom_id}{"potential_energy"}
    	&& $atom_site->{$atom_id}{"potential_energy"} > $cutoff;
    }

    return $atom_site;
}

# ------------------------- Various potential functions ----------------------- #

#
# Hard sphere potential function. Described as:
#     0,   r_{ij} >= vdw_{i} + vdw_{j}
#     Inf, r_{ij} <  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            Inf - infinity;
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     two values: 0 or "Inf" (infinity).
#

sub hard_sphere
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance < $vdw_length ** 2 ) {
	return "Inf";
    } else {
	return 0;
    }
}

#
# Soft sphere potential function. Described as:
#     epsilon * ( vdw_{i} + vdw_{j} / r_{ij} ) ** n , r_{ij} <= vdw_{i} + vdw_{j}
#     0, r_{ij} >  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            epsilon - energy coefficient; TODO: should study more about it.
#            n - number increases the slope of potential.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by soft sphere potential.
#

sub soft_sphere
{
    my ( $target_atom, $neighbour_atom, $epsilon, $n ) = @_;

    $epsilon //= 1.0;
    $n //= 12;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance <= $vdw_length ** 2 ) {
	return $epsilon * ( $vdw_length / sqrt( $distance ) )**$n;
    } else {
	return 0;
    }
}


#
# Exponential potential function. Described as:
#     epsilon * exp( - ( r_{ij} / vdw_{i} + vdw_{j} ) ) ** m ,
#        r_{ij} <= vdw_{i} + vdw_{j}
#     0, r_{ij} >  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            epsilon - energy coefficient; TODO: should study more about it.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by exponential potential.
#

sub exponential
{
    my ( $target_atom, $neighbour_atom, $epsilon, $m ) = @_;

    $epsilon //= 1.0;
    $m //= 1.0;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance <= $vdw_length ** 2 ) {
	return $epsilon * exp( - ( sqrt( $distance ) / $vdw_length ) ** $m );
    } else {
	return 0;
    }
}

#
# Leonard-Jones potential function. Described as:
#
# 4 * epsilon * [ ( sigma / r ) ** 12 - ( sigma / r ) ** 6 ]
#
#     where: r - distance between center of atoms;
#            epsilon - energy coefficient;  TODO: should study more about it.
#            sigma - sum of Van der Waals radii of two atoms.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by Leonard-Jones potential.
#

sub leonard_jones
{
    my ( $target_atom, $neighbour_atom, $epsilon ) = @_;

    $epsilon //= -1.0;

    my $sigma =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $r =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $r <= ( $sigma ** 2.5 ) ** 2 ) { # TODO:should smoothen function cutoff.
	return 4 * $epsilon * ( ( $sigma / $r ) ** 12 - ( $sigma / $r ) ** 6 );
    } else {
	return 0;
    }
}

1;
