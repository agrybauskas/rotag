package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( max );

use lib qw( ./ );
use ConnectAtoms qw( grid_box );
use LoadParams qw( covalent_radii
                   vdw_radii );
use Data::Dumper;

my $covalent_file = "../../parameters/covalent_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Parameters.
#

my %COVALENT_RADII = %{ covalent_radii( $covalent_file ) };

my $MAX_BOND_LENGTH =
    max( map { @{ $COVALENT_RADII{$_}{"bond_length"} } }
	 keys( %COVALENT_RADII ) ) * 2;

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

    # Identifies atoms that are in a clash with other atoms.

    # For each cell, checks neighbouring cells.
    my %atom_clashes = %{ $atom_site };
    my @cell_idx;

    # Creates box around atoms, makes grid with edge length of max covalent radii.
    my $grid_box =
    	grid_box( $atom_site, $MAX_BOND_LENGTH );

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

    #     # Atoms that have been already checked for connections.
    # 	my @checked_atoms;

    # 	# Checks, if there are clashes between atoms.
    # 	foreach my $atom_id ( @{ $grid_box->{$cell} } ) {
    # 	    push( @checked_atoms, $atom_id ); # Marks as visited atom.
    # 	    foreach my $neighbour_id ( @neighbour_cells ) {
    # 		if( check_distance($atom_site->{"data"}{"$atom_id"},
    # 				   $atom_site->{"data"}{"$neighbour_id"} )
    # 		    eq "clash" ) {
    # 		    push( @{ $atom_clashes{"data"}
    # 			                  {$atom_id}
    # 			                  {"connections"} },
    # 		          $neighbour_id );
    # 		}
    # 	    }
    # 	}
    # }

    # return \%atom_clashes;

}

1;
