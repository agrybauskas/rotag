package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_only );

use lib qw( ./ );
use CifParser qw( select_atom_data filter_atoms );
use AlterMolecule qw( bond_torsion );
use LinearAlgebra qw( evaluate_matrix
                      matrix_product
                      switch_ref_frame
                      vectorize );
use LoadParams qw( rotatable_bonds );

my $parameter_file = "../../parameters/rotatable_bonds.csv";

# ------------------------ Idealistic sidechain models ------------------------ #

#
# Parameters.
#

my %ROTATABLE_BONDS = %{ rotatable_bonds( $parameter_file ) };

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final. One model per
# residue.
#

#
# Model that uses only rotation around single bonds.
# Input  (1 arg):
# Output (1 arg):
#

sub rotation_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Selects specified atom(s) id(s).
    my $target_atom_id =
	&select_atom_data( [ "id" ],
			   &filter_atoms( $atom_specifier,
					  $atom_site ) );

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

    for my $id ( @$target_atom_id ) {
    	$residue_name = $atom_site->{"data"}{"@$id"}{"label_comp_id"};
    	$atom_type = $atom_site->{"data"}{"@$id"}{"label_atom_id"};
    	@rotatable_bonds = @{ $ROTATABLE_BONDS{$residue_name}{$atom_type} };
    	@transf_matrices = ();
    	$target_atom_coord =
    	    &select_atom_data( [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    			       &filter_atoms( $atom_specifier,
    					      $atom_site ) );

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
    		    &select_atom_data(
    		    [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		    &filter_atoms(
    			{ "label_atom_id" => [ "N" ] },
    			$atom_site ) );
    	    } else {
    		$side_atom_coord =
    		    &select_atom_data(
    		    [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		    &filter_atoms(
    			{ "label_atom_id" => [ $rotatable_bonds[$i-1][1] ] },
    			$atom_site ) );
    	    }

    	    $mid_atom_coord =
    		&select_atom_data(
    		[ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		&filter_atoms(
    		    { "label_atom_id" => [ $mid_atom_type ] },
    		    $atom_site ) );

    	    $up_atom_coord =
    		&select_atom_data(
    		[ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		&filter_atoms(
    		    { "label_atom_id" => [ $up_atom_type ] },
    		    $atom_site ) );

    	    # Creates and appends matrices to a list of matrices that later
    	    # will be multiplied.
    	    push( @transf_matrices,
    		  bond_torsion( $angle_symbol,
				@$mid_atom_coord,
    				@$up_atom_coord,
    				@$side_atom_coord ) );
    	}

    	$atom_site->{"data"}{"@$id"}{"conformation"} =
    	    matrix_product( @transf_matrices,
    			    vectorize( @{ $target_atom_coord } ) );
    }

    return $atom_site;
}

1;
