package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_only );

use lib qw( ./ );
use AlterMolecule qw( bond_torsion );
use LinearAlgebra qw( evaluate_matrix
                      matrix_product
                      switch_ref_frame
                      vectorize );
use MoleculeProperties qw( %ROTATABLE_BONDS );
use PDBxParser qw( filter_atoms select_atom_data );

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

    # iterates through target atom(s) and assigns conformational equations which
    # can produce pseudo-atoms later.
    for my $id ( keys %{ $atom_site } ) {
	my $residue_id =   $atom_site->{"$id"}{"label_seq_id"};
    	my $residue_name = $atom_site->{"$id"}{"label_comp_id"};
    	my $atom_type =    $atom_site->{"$id"}{"label_atom_id"};
	my @atom_coord = ( $atom_site->{"$id"}{"Cartn_x"},
			   $atom_site->{"$id"}{"Cartn_y"},
			   $atom_site->{"$id"}{"Cartn_z"} );

	# Checks, if according to residue name and atom type, conformational
	# model can be applied.
	my @rotatable_bonds;
	if( exists $ROTATABLE_BONDS{$residue_name}{$atom_type} ) {
	    @rotatable_bonds =
		@{ $ROTATABLE_BONDS{$residue_name}{$atom_type} };
	} else {
	    next;
	}

    	my @transf_matrices; # Matrices for transforming atom coordinates.

    	for( my $i = 0; $i < scalar( @rotatable_bonds ) - 1; $i++ ) {
    	    my $mid_atom_type = $rotatable_bonds[$i][0];
    	    my $up_atom_type = $rotatable_bonds[$i][1];
	    $mid_atom_type =~ s/\s//g;
    	    $up_atom_type =~ s/\s//g;

    	    my $angle_symbol = "chi${i}";

    	    # Information about side atom is stored in rotatable bonds array,
    	    # except for CA atom.
	    my $side_atom_coord;
    	    if( $mid_atom_type eq "CA" ) {
    	    	$side_atom_coord =
    	    	    select_atom_data(
    	    	    filter_atoms( $atom_site,
    	    	    { "label_seq_id"  => [ "$residue_id" ],
		      "label_atom_id" => [ "N" ] } ),
    	    	    [ "Cartn_x", "Cartn_y", "Cartn_z" ] );
    	    } else {
    		$side_atom_coord =
    		    select_atom_data(
    		    filter_atoms( $atom_site,
    		    { "label_atom_id" => [ $rotatable_bonds[$i-1][1] ] } ),
    		    [ "Cartn_x", "Cartn_y", "Cartn_z" ] );
    	    }

    	    my $mid_atom_coord =
    		select_atom_data(
    		filter_atoms( $atom_site,
    		{ "label_atom_id" => [ $mid_atom_type ],
		  "label_seq_id"  => [ "$residue_id" ] } ),
    		[ "Cartn_x", "Cartn_y", "Cartn_z" ] );

    	    my $up_atom_coord =
    		select_atom_data(
    		filter_atoms( $atom_site,
    		{ "label_atom_id" => [ $up_atom_type ],
		  "label_seq_id"  => [ "$residue_id" ] } ),
    		[ "Cartn_x", "Cartn_y", "Cartn_z" ] );

    	    # Creates and appends matrices to a list of matrices that later
    	    # will be multiplied.
    	    push( @transf_matrices,
    		  bond_torsion( @{ $mid_atom_coord },
    				@{ $up_atom_coord },
    				@{ $side_atom_coord },
    				$angle_symbol ) );
    	}

    	$atom_site->{"$id"}{"conformation"} =
    	    matrix_product( @transf_matrices,
    			    vectorize( \@atom_coord ) );
    }

    return $atom_site;
}

1;
