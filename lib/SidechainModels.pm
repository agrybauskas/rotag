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
use PDBxParser qw( filter );

# ------------------------ Idealistic sidechain models ------------------------ #

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
	uniq( @{ filter( { "atom_site" => $atom_site,
			   "data" => [ "label_seq_id" ],
			   "is_list" => 1 } ) } );

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_id ( @residue_ids ) {
    	my $residue_site =
	    filter( { "atom_site" => $atom_site,
		      "include" => { "label_seq_id" => [ $residue_id ] } } );
	my $side_chain_site =
	    filter( { "atom_site" => $residue_site,
		      "exclude" => { "label_atom_id" =>
				       [ "N", "C", "O", "H", "H2", "HA" ] } } ) ;

    	# Determines sidechain atoms, CA atom and then finds rotatable bonds.
    	connect_atoms( $residue_site );
    	hybridization( $residue_site );

    	# Determines rotatable bonds.
	my $ca_id =
	    filter( { "atom_site" => $side_chain_site,
		      "label_atom_id" => [ "CA" ],
		      "data" => [ "id" ],
		      "is_list" => 1 } );
	my $cb_id =
	    filter( { "atom_site" => $side_chain_site,
		      "label_atom_id" => [ "CB" ],
		      "data" => [ "id" ],
		      "is_list" => 1 } );
    	my $rotatable_bonds =
	    rotatable_bonds( $side_chain_site, @{ $ca_id }, @{ $cb_id } );

	# Assigns names to rotatable bond angles.
	# TODO: check if later names of rotation_only dihedral angles match
	# all_dihedral values.
	my @rotatable_bonds;
	my %rot_bond_names;
	my $angle_id = 0;
	for my $atom_id (
	    sort { scalar @{ $rotatable_bonds->{$a} }
	       <=> scalar @{ $rotatable_bonds->{$a} } }
	         keys %{ $rotatable_bonds } ) {
    	    # Filters out redundant rotatable bonds.
    	    for my $rotatable_bond ( @{ $rotatable_bonds->{$atom_id} } ) {
    		if( ! grep { $rotatable_bond->[0] eq $_->[0]
    			  && $rotatable_bond->[1] eq $_->[1]} @rotatable_bonds ){
    		    push( @rotatable_bonds, $rotatable_bond );
		    $rot_bond_names{"$rotatable_bond->[0],$rotatable_bond->[1]"}=
			"chi${angle_id}";
		    $angle_id++;
    		}
	    }
	}

    	for my $atom_id ( keys %{ $residue_site }  ) {
    	    my @atom_coord = ( $atom_site->{"$atom_id"}{"Cartn_x"},
    			       $atom_site->{"$atom_id"}{"Cartn_y"},
    			       $atom_site->{"$atom_id"}{"Cartn_z"} );

    	    if( ! exists $rotatable_bonds->{$atom_id} ) { next; }

    	    my @transf_matrices; # Matrices for transforming atom coordinates.

    	    for my $rotatable_bond ( @{ $rotatable_bonds->{$atom_id} } ) {
		# First, checks if rotatable bond has fourth atom produce
		# dihedral angle. It is done by looking at atom connections - if
		# rotatable bond ends with terminal atom, then this bond is
		# excluded.
		if( scalar( @{ $residue_site->{$rotatable_bond->[1]}
			       {"connections"} } ) < 2 ){ next; }

    		my $mid_atom_id = $rotatable_bond->[0];
    		my $up_atom_id  = $rotatable_bond->[1];
    		my @mid_connections = # Excludes up atom.
    		    grep { $_ ne $up_atom_id }
    		    @{ $residue_site->{$mid_atom_id}{"connections"} };
    		my @mid_connection_names = # Excludes up atom.
    		    map { $residue_site->{$_}{"label_atom_id"} }
    		    @mid_connections;
    		my $side_atom_name =
    		    sort_by_priority( \@mid_connection_names )->[0];
    		my $side_atom_id =
		    filter( { "atom_site" => $residue_site,
			      "include" =>
			    { "label_atom_id" => [ $side_atom_name ] },
			      "data" => [ "id" ],
			      "is_list" => 1 } )->[0];

		my $angle_symbol;
		if( exists $rot_bond_names{"$mid_atom_id,$up_atom_id"} ) {
		    $angle_symbol =
			$rot_bond_names{"$mid_atom_id,$up_atom_id"};
		} elsif( exists $rot_bond_names{"$up_atom_id,$mid_atom_id"} ) {
		    $angle_symbol =
			$rot_bond_names{"$up_atom_id,$mid_atom_id"};
		} else {
		    die( "No corresponding angle names were found" );
		}

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
