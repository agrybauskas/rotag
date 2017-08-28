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
use LoadParams qw( rotatable_bonds );
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
#     $atom_specifier - data structure that specifies atom (see PDBxParser.pm)
# Output:
#     $atom_site - modified $atom_site with added equation describing
#     conformational space.
#

sub rotation_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Selects specified atom(s) id(s).
    my $target_atom_id =
	select_atom_data( filter_atoms( $atom_site, $atom_specifier ),
			  [ "id" ] );

    # Iterates through target atom(s) and assigns conformational equations which
    # can produce pseudo-atoms later.
    my $residue_name;
    my $atom_type;

    my @rotatable_bonds;

    my $mid_atom_type;
    my $up_atom_type;

    my $mid_atom_coord;
    my $up_atom_coord;
    my $side_atom_coord;
    my $target_atom_coord;

    my @transf_matrices; # Matrices for transforming atom coordinates.

    for my $id ( @{ $target_atom_id } ) {
    	$residue_name = $atom_site->{"@$id"}{"label_comp_id"};
    	$atom_type = $atom_site->{"@$id"}{"label_atom_id"};
    	@rotatable_bonds = @{ rotatable_bonds->{$residue_name}{$atom_type} };
    	undef @transf_matrices;
    	$target_atom_coord =
    	    select_atom_data( filter_atoms( $atom_site, $atom_specifier ),
			      [ "Cartn_x", "Cartn_y", "Cartn_z" ] );


    	# Creates matrices for atom alterations.
    	my $angle_symbol; # Because side chain might have multiple
    	                  # rotatable bonds, there must be distinct
    	                  # symbols for different dihedral angles.

    	for( my $i = 0; $i < scalar( @rotatable_bonds ) - 1; $i++ ) {
    	    $angle_symbol = "chi${i}";
    	    $mid_atom_type = $rotatable_bonds[$i][0];
    	    $mid_atom_type =~ s/\s//g;
    	    $up_atom_type = $rotatable_bonds[$i][1];
    	    $up_atom_type =~ s/\s//g;

    	    # Information about side atom is stored in rotatable bonds array,
    	    # except for CA atom.
    	    if( $mid_atom_type eq "CA" ) {
    		$side_atom_coord =
    		    select_atom_data(
		    filter_atoms( $atom_site,
		    { "label_atom_id" => [ "N" ] } ),
		    [ "Cartn_x", "Cartn_y", "Cartn_z" ] );

    	    } else {
    		$side_atom_coord =
    		    select_atom_data(
		    filter_atoms( $atom_site,
		    { "label_atom_id" => [ $rotatable_bonds[$i-1][1] ] } ),
		    [ "Cartn_x", "Cartn_y", "Cartn_z" ] );
    	    }

    	    $mid_atom_coord =
		select_atom_data(
		filter_atoms( $atom_site,
		{ "label_atom_id" => [ $mid_atom_type ] } ),
		[ "Cartn_x", "Cartn_y", "Cartn_z" ] );

    	    $up_atom_coord =
		select_atom_data(
		filter_atoms( $atom_site,
		{ "label_atom_id" => [ $up_atom_type ] } ),
		[ "Cartn_x", "Cartn_y", "Cartn_z" ] );

    	    # Creates and appends matrices to a list of matrices that later
    	    # will be multiplied.
    	    push( @transf_matrices,
    		  bond_torsion( @{ $mid_atom_coord },
    				@{ $up_atom_coord },
    				@{ $side_atom_coord },
				$angle_symbol ) );
    	}

    	$atom_site->{"@$id"}{"conformation"} =
    	    matrix_product( @transf_matrices,
    			    vectorize( @{ $target_atom_coord } ) );
    }

    return $atom_site;
}

1;
