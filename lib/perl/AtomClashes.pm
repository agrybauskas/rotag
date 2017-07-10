package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( max );

use lib qw( ./ );
use CifParser qw( filter_atoms
                  select_atom_data );
use ConnectAtoms qw( check_distance
                     grid_box );
use LoadParams qw( vdw_radii );
use Data::Dumper;
my $vdw_file = "../../parameters/vdw_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Parameters.
#

my %VDW_RADII = %{ vdw_radii( $vdw_file ) };

my $MAX_BOND_LENGTH =
    max( map { $VDW_RADII{$_} }
	 keys( %VDW_RADII ) ) * 2;

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

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

    # Identifies atoms that are in a clash with other atoms.

    # For each cell, checks neighbouring cells.
    my %atom_clashes = %{ $atom_site };
    my @cell_idx;

    # Creates box around atoms, makes grid with edge length of max covalent radii.
    my $grid_box =
    	grid_box( $atom_site, $MAX_BOND_LENGTH );

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

    	# Checks, if there are clashes between atoms.
    	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
	    if( $atom_id ~~ @spec_atom_ids ) {
		foreach my $neighbour_id ( @neighbour_cells ) {
		    if( not $neighbour_id ~~ @spec_atom_ids ) {
			# TODO: make check_distance dependent on parameters.
			if( check_distance(
				$atom_site->{"data"}{"$atom_id"},
				$atom_site->{"data"}{"$neighbour_id"} ) eq "clash"
			    && $atom_id != $neighbour_id ) {
			    push( @{ $atom_clashes{"data"}
				     {$atom_id}
				     {"clashes"} },
				  $neighbour_id );
			}
		    }
		}
	    }
    	}
    }

    return \%atom_clashes;
}

1;
