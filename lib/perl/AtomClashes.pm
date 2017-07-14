package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( any max );

use lib qw( ./ );
use CifParser qw( filter_atoms select_atom_data );
use ConnectAtoms qw( connect_atoms grid_box );
use LoadParams qw( rotatable_bonds vdw_radii );
use Data::Dumper;

my $rot_bonds_file = "../../parameters/rotatable_bonds.csv";
my $vdw_file = "../../parameters/vdw_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

my %ROTATABLE_BONDS = %{ rotatable_bonds( $rot_bonds_file ) };

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

    if( $distance < $vdw_length ** 2 ) {
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

    # Identifies atoms that are in a clash with other atoms.

    # For each cell, checks neighbouring cells.
    my %atom_clashes = %{ connect_atoms( $atom_site ) };
    my @cell_idx;

    # Creates box around atoms, makes grid with edge length of max covalent radii.
    my $grid_box =
    	grid_box( $atom_site, $MAX_VDW_RADIUS );

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
    	    if( any { $atom_id eq $_ } @spec_atom_ids ) {
    		foreach my $neighbour_id ( @neighbour_cells ) {
    		    if( not any { $neighbour_id eq $_ } @spec_atom_ids ) {
    			if( is_colliding( $atom_site->{"data"}{"$atom_id"},
    					  $atom_site->{"data"}{"$neighbour_id"} )
    			    && $atom_id ne $neighbour_id ) {
    			    push( @{ $atom_clashes{"data"}
    				     {$atom_id}
    				     {"clashes"} },
    				  $neighbour_id );
    			}
    		    }
    		}
    	    }
    	}

	# Removes clashes, if there are more than two.
	my $res_id;
	my $clash_res_id;
	my $atom_name;
	my $clash_name;
	my $second_neighbour;

	foreach my $atom_id ( keys %{ $atom_clashes{"data"} } ) {
	    $res_id = $atom_clashes{"data"}{$atom_id}{"label_seq_id"};
	    $atom_name = $atom_clashes{"data"}{$atom_id}{"label_atom_id"};
	    foreach my $clash_id (
		@{ $atom_clashes{"data"}{$atom_id}{"clashes"} } ) {
		$clash_res_id = $atom_clashes{"data"}{$clash_id}{"label_seq_id"};
		$clash_name = $atom_clashes{"data"}{$clash_id}{"label_atom_id"};
		if( $clash_res_id eq $res_id ) {
		    # Checks if clashing atom is not second after atom that
		    # currently is connected to.
		    for my $i (
			@{ $atom_clashes{"data"}{$atom_id}{"connections"} } ) {
		    for my $j (
			@{ $atom_clashes{"data"}{$i}{"connections"} } ) {
		        if( $clash_name eq
			    $atom_clashes{"data"}{$j}{"label_atom_id"} ) {
			    last;
			}
		    } last }
		    # $clash_name = $atom_clashes{"data"}{$clash_id}{"label_seq_id"};
		}
	    }
	}
    # 	foreach my $atom_id ( keys %{ $atom_clashes{"data"} } ) {
    # 	    if( exists $atom_clashes{"data"}{$atom_id}{"clashes"}
    # 	     && scalar( @{ $atom_clashes{"data"}{$atom_id}{"clashes"} } ) > 3 ) {
    # 		delete $atom_clashes{"data"}{$atom_id};
    # 	    }
    # 	}
    }

    # return \%atom_clashes;
}

1;
