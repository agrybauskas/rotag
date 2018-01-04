package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_only );

use List::MoreUtils qw( uniq );

use AlterMolecule qw( bond_torsion );
use ConnectAtoms qw( connect_atoms
                     hybridization
                     rotatable_bonds
                     sort_by_priority );
use LinearAlgebra qw( flatten
                      mult_matrix_product
                      reshape );
use PDBxParser qw( filter_atoms
                   select_atom_data );

# ------------------------ Idealistic sidechain models ------------------------ #

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final. One model per
# residue.
#

#
# Model that uses only rotation around single bonds.
# Input:
#     $atom_site - atom data structure.
# Output:
#     $atom_site - modified $atom_site with added equation describing
#     conformational space.
#

sub rotation_only
{
    my ( $atom_site ) = @_;

    # Determines all residue ids present in atom site.
    my @residue_ids =
	uniq( @{ flatten( [ select_atom_data( $atom_site,
					      [ "label_seq_id" ] ) ] ) } );

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_id ( @residue_ids ) {
	my $residue_site =
	    filter_atoms( $atom_site, { "label_seq_id" => [ $residue_id ] } );
	my @atom_ids = keys %{ $residue_site };

	# Determines sidechain atoms, CA atom and then finds rotatable bonds.
	connect_atoms( $residue_site );
	hybridization( $residue_site );

	my @main_chain_ids = # Main chain atoms except for CA.
	    @{ flatten(
	       [ select_atom_data(
	         filter_atoms(
	         $residue_site,
	         { "label_atom_id" => [ "N", "C", "O", "OXT", "H", "H2", "H3",
					"HXT" ] } ),
	         [ "id" ] ) ] ) };
	my @side_chain_ids;
	for my $atom_id ( @atom_ids ) {
	    if( ! grep { $atom_id eq $_ } @main_chain_ids ) {
		push( @side_chain_ids, $atom_id );
	    }
	}
	my $ca_id =
	    select_atom_data(
	    filter_atoms(
	    $residue_site,
	    { "label_atom_id" => [ "CA" ] } ),
	    [ "id" ] )->[0][0];

	my $rotatable_bonds =
	    rotatable_bonds(
		filter_atoms( $residue_site,
			      { "id" => \@side_chain_ids } ),
		$ca_id );

	for my $atom_id ( @atom_ids ) {
	    my @atom_coord = ( $atom_site->{"$atom_id"}{"Cartn_x"},
			       $atom_site->{"$atom_id"}{"Cartn_y"},
			       $atom_site->{"$atom_id"}{"Cartn_z"} );

	    if( ! exists $rotatable_bonds->{$atom_id} ) { next; }

	    # Determines have one rotational bond and, also, part of it. Those
	    # atoms are excluded from analysis.
	    if( ( scalar(  @{ $rotatable_bonds->{$atom_id} } ) == 1 )
	     && ( $atom_id == $rotatable_bonds->{$atom_id}->[0][0]
	       || $atom_id == $rotatable_bonds->{$atom_id}->[0][1] ) ) { next; }

	    my @rotatable_bonds = @{ $rotatable_bonds->{$atom_id} };

	    my @transf_matrices; # Matrices for transforming atom coordinates.

	    for( my $i = 0; $i < scalar( @rotatable_bonds ) - 1; $i++ ) {
		my $mid_atom_id = $rotatable_bonds[$i][0];
		my $up_atom_id  = $rotatable_bonds[$i][1];

		my @mid_connections = # Excludes up atom.
		    grep { $_ ne $up_atom_id }
		    @{ $residue_site->{$mid_atom_id}{"connections"} };
		my @mid_connection_names = # Excludes up atom.
		    map { $residue_site->{$_}{"label_atom_id"} }
		    @mid_connections;
		my $side_atom_name =
		    sort_by_priority( \@mid_connection_names )->[0];
		my $side_atom_id =
		    select_atom_data(
		    filter_atoms(
		    $residue_site,
		    { "label_atom_id" => [ $side_atom_name ] } ),
		    [ "id" ] )->[0][0];

		# TODO: check if angle names do not change due to different
		# ordering of atom ids.
		my $angle_symbol = "chi${i}";

		my $mid_atom_coord =
		    [ $residue_site->{$mid_atom_id}{"Cartn_x"},
		      $residue_site->{$mid_atom_id}{"Cartn_y"},
		      $residue_site->{$mid_atom_id}{"Cartn_z"} ];
		my $up_atom_coord =
		    [ $residue_site->{$up_atom_id}{"Cartn_x"},
		      $residue_site->{$up_atom_id}{"Cartn_y"},
		      $residue_site->{$up_atom_id}{"Cartn_z"} ];
		my $side_atom_coord =
		    [ $residue_site->{$side_atom_id}{"Cartn_x"},
		      $residue_site->{$side_atom_id}{"Cartn_y"},
		      $residue_site->{$side_atom_id}{"Cartn_z"} ];

    	    # Creates and appends matrices to a list of matrices that later
    	    # will be multiplied.
    	    push( @transf_matrices,
    	    	  @{ bond_torsion( $mid_atom_coord,
				   $up_atom_coord,
				   $side_atom_coord,
				   $angle_symbol ) } );
	    }

	    $atom_site->{$atom_id}{"conformation"} =
		mult_matrix_product(
		    [ @transf_matrices,
		      @{ reshape( [ @atom_coord, 1 ], [ 4, 1 ] ) } ] );
	}
    }

    return $atom_site;
}

1;
