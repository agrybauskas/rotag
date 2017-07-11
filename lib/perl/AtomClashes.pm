package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( max );

use lib qw( ./ );
use CifParser qw( filter_atoms
                  select_atom_data );
use ConnectAtoms qw( is_connected
                     grid_box );
use LoadParams qw( vdw_radii );
use Data::Dumper;
my $vdw_file = "../../parameters/vdw_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

my %VDW_RADII = %{ vdw_radii( $vdw_file ) };

my $MAX_VDW_RADIUS =
    max( map { $VDW_RADII{$_} }
	 keys( %VDW_RADII ) ) * 2;

sub is_colliding
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $vdw_length =
	$VDW_RADII{$target_atom->{"type_symbol"}}
      + $VDW_RADII{$neighbour_atom->{"type_symbol"}};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    # Checks, if distance between atom pairs is in one of the combinations.
    my $is_colliding;

    if( ! ( is_connected( $target_atom, $neighbour_atom ) )
     && ! ( $distance > $vdw_length ** 2 ) ) {
	$is_colliding = 1;
    } else {
	$is_colliding = 0;
    }

    return $is_colliding;
}

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
#

sub radius_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Clashes of all atoms analyzed, if no specific atoms are selected.
    $atom_specifier = { "group_pdb" => [ "ATOM" ] } unless $atom_specifier;
    my @spec_atom_ids = # Atom ids selected by $atom_specifier.
    	map { $_->[0] }
        @{ select_atom_data( [ "id" ],
    			     filter_atoms( $atom_specifier, $atom_site ) ) };

    # # Identifies atoms that are in a clash with other atoms.

    # # For each cell, checks neighbouring cells.
    # my %atom_clashes = %{ $atom_site };
    # my @cell_idx;

    # # Creates box around atoms, makes grid with edge length of max covalent radii.
    # my $grid_box =
    # 	grid_box( $atom_site, $MAX_VDW_RADIUS );

    # # Checks for neighbouring cells for each cell.
    # foreach my $cell ( keys %{ $grid_box } ) {
    # 	@cell_idx = split( ",", $cell );
    # 	my @neighbour_cells; # The array will contain all atoms of the
    #                          # neighbouring 26 cells.

    # 	# $i represents x, $j - y, $k - z coordinates.
    # 	for my $i ( ( $cell_idx[0] - 1..$cell_idx[0] + 1 ) ) {
    # 	for my $j ( ( $cell_idx[1] - 1..$cell_idx[1] + 1 ) ) {
    # 	for my $k ( ( $cell_idx[2] - 1..$cell_idx[2] + 1 ) ) {
    # 	if( exists $grid_box->{"$i,$j,$k"} ) {
    # 	    push( @neighbour_cells, @{ $grid_box->{"$i,$j,$k"} } ); } } } }

    # 	# Checks, if there are clashes between atoms.
    # 	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
    # 	    if( ! grep { $atom_id } @spec_atom_ids ) {
    # 		foreach my $neighbour_id ( @neighbour_cells ) {
    # 		    if( ! grep { $neighbour_id } @spec_atom_ids ) {
    # 			if( is_colliding( $atom_site->{"data"}{"$atom_id"},
    # 					  $atom_site->{"data"}{"$neighbour_id"} )
    # 			    && $atom_id != $neighbour_id ) {
    # 			    push( @{ $atom_clashes{"data"}
    # 				     {$atom_id}
    # 				     {"clashes"} },
    # 				  $neighbour_id );
    # 			}
    # 		    }
    # 		}
    # 	    }
    # 	}
    # }

    # return \%atom_clashes;
}

1;
