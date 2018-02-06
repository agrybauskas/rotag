package Measure;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( all_dihedral
                     bond_angle
                     bond_length
                     dihedral_angle
                     rmsd );

use Math::Trig;
use List::MoreUtils qw( uniq );

use ConnectAtoms qw( connect_atoms
                     rotatable_bonds
                     sort_by_priority );
use PDBxParser qw( filter );
use LinearAlgebra qw( matrix_sub
                      vector_cross );
use Data::Dumper;
# ----------------------------- Molecule parameters --------------------------- #

#
# Calculates various parameters that describe molecule or atoms, such as, bond
# length, dihedral angle, torsion angle, RMSD and etc.
#

#
# Calculates bond length of given two atoms.
# Input:
#     $atom_coord - two matrices of x, y, z coordinates corresponding  two atoms.
# Output:
#     $bond_length - length of the bond (in angstroms).
#

sub bond_length
{
    my ( $atom_coord ) = @_;

    my $bond_length =
	sqrt( ( $atom_coord->[1][0] - $atom_coord->[0][0] )**2
	    + ( $atom_coord->[1][1] - $atom_coord->[0][1] )**2
	    + ( $atom_coord->[1][2] - $atom_coord->[0][2] )**2 );

    return $bond_length;
}

#
# Calculates angle between three atoms.
# Input:
#     $atom_coord - three matrices of x,y,z coordinates corresponding to three
#     atoms.
# Output:
#     $bond_angle - angle in radians.
#

sub bond_angle
{
    my ( $atom_coord ) = @_;

    # Angle between three atoms (in radians) in 3-D space can be calculated by
    # the formula:
    #                            ->   ->      ->         ->
    #            theta = arccos( AB * BC / || AB || * || BC || )

    # This formula is applied to atoms where vectors are the substraction of
    # coordinates of two atoms. Suppose, one of the side atom is A, B - middle
    # and C - remaining atom. Order of side atoms is irrelevant.
    my @vector_ab = ( $atom_coord->[0][0] - $atom_coord->[1][0],
		      $atom_coord->[0][1] - $atom_coord->[1][1],
		      $atom_coord->[0][2] - $atom_coord->[1][2] );
    my @vector_bc = ( $atom_coord->[2][0] - $atom_coord->[1][0],
		      $atom_coord->[2][1] - $atom_coord->[1][1],
		      $atom_coord->[2][2] - $atom_coord->[1][2] );

    my $vector_product = $vector_ab[0] * $vector_bc[0] +
	               + $vector_ab[1] * $vector_bc[1]
	               + $vector_ab[2] * $vector_bc[2];

    my $length_ab = sqrt( $vector_ab[0]**2
			+ $vector_ab[1]**2
			+ $vector_ab[2]**2 );
    my $length_bc = sqrt( $vector_bc[0]**2
			+ $vector_bc[1]**2
			+ $vector_bc[2]**2 );

    my $bond_angle = acos( $vector_product / ( $length_ab * $length_bc ) );

    return $bond_angle;
}

#
# Calculates dihedral angle of four given atoms.
# Input:
#     $atom_coord - four matrices of x,y,z coordinates corresponding to four
#     atoms.
# Output:
#     $dihedral - dihedral angle in radians.
#

sub dihedral_angle
{
    my ( $atom_coord ) = @_;

    #                  -> ->    ->
    # Creates vectors: a, b and c, that are translated to global reference frame.
    # Picture of vectors:
    #                                   ->  O ->
    #                                   b  /  c
    #                              -> CA---C
    #                              a /
    #                               N
    my $vector_a = matrix_sub( [ $atom_coord->[1] ], [ $atom_coord->[0] ] );
    my $vector_b = matrix_sub( [ $atom_coord->[2] ], [ $atom_coord->[1] ] );
    my $vector_c = matrix_sub( [ $atom_coord->[3] ], [ $atom_coord->[2] ] );

    #                                               ->    -> ->    ->
    # Creates normal vectors from the vector pairs: a and b, b and c.
    my $vector_cross_ab = vector_cross( @$vector_a, @$vector_b );
    my $vector_cross_bc = vector_cross( @$vector_b, @$vector_c );

    # Calculates length for each cross product.
    my $vector_length_ab = sqrt( $vector_cross_ab->[0]**2
			       + $vector_cross_ab->[1]**2
			       + $vector_cross_ab->[2]**2 );
    my $vector_length_bc = sqrt( $vector_cross_bc->[0]**2
			       + $vector_cross_bc->[1]**2
			       + $vector_cross_bc->[2]**2 );

    # Calculates normal vectors for each cross product of two vectors.
    my @normal_vector_ab = map { $_ / $vector_length_ab } @{ $vector_cross_ab };
    my @normal_vector_bc = map { $_ / $vector_length_bc } @{ $vector_cross_bc };

    # Finishes orthonormal frame from normal vector ab, vector b and its cross
    # product.
    my $vector_length_b = sqrt( $vector_b->[0][0]**2
    			      + $vector_b->[0][1]**2
    			      + $vector_b->[0][2]**2 );
    my @normal_vector_b = map { $_ / $vector_length_b } @{ $vector_b->[0] };

    my @orthonormal_cross =
	vector_cross( \@normal_vector_ab, \@normal_vector_b );

    # Using orthonormal frame, projections from vector a and c are
    # generated and angle calculated.
    my $dihedral_angle =
	- atan2( $orthonormal_cross[0][0] * $normal_vector_bc[0]
	       + $orthonormal_cross[0][1] * $normal_vector_bc[1]
	       + $orthonormal_cross[0][2] * $normal_vector_bc[2],
	         $normal_vector_ab[0] * $normal_vector_bc[0]
	       + $normal_vector_ab[1] * $normal_vector_bc[1]
	       + $normal_vector_ab[2] * $normal_vector_bc[2] );

    return $dihedral_angle;
}

#
# Calculates dihedral angles for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it).
# Usage of connect_atoms and hybridization functions are necessary for correct
# calculations.
# Input:
#     $atom_site - atom data structure.
# Output:
#     $residue_angles - data structure that relates residue id and angle values.
#     Ex.:
#       resi_id, angle_name, angle value
#     { 18 => { 'chi0' => '-3.14' } }
#

sub all_dihedral
{
    my ( $atom_site ) = @_;

    # Collects non-redundant ids of given amino acid residues.
    my @residue_ids =
	uniq( @{ filter( { "atom_site" => $atom_site,
			   "data" => [ "label_seq_id" ],
			   "is_list" => 1 } ) } );

    # TODO: try to figure out, where these functions that modify states should be
    # used during the run of functions (in the beginning of script or during
    # every execute of the function that needs the result of the previously
    # mentioned functions).
    my %atom_site = %{ $atom_site }; # Copy of $atom_site.

    connect_atoms( \%atom_site );

    # Iterates through residue ids and, according to the parameter file,
    # calculates dihedral angles of each rotatable bond.
    my %residue_angles;

    for my $residue_id ( @residue_ids ) {
    	my $residue_site =
    	    filter( { "atom_site" => \%atom_site,
    		      "include" => { "label_seq_id" => [ $residue_id ] } } );

    	my $rotatable_bonds = rotatable_bonds( $residue_site );
	my %uniq_rotatable_bonds; # Unique rotatable bonds.
	for my $atom_id ( keys %{ $rotatable_bonds } ) {
	    for my $angle_name ( keys %{ $rotatable_bonds->{"$atom_id"} } ){
		if( ! exists $uniq_rotatable_bonds{"$angle_name"} ) {
		    $uniq_rotatable_bonds{"$angle_name"} =
			$rotatable_bonds->{"$atom_id"}{"$angle_name"};
		}
	    }
	}

    	# Calculates every dihedral angle.
    	my %angle_values;

    	for my $angle_name ( keys %uniq_rotatable_bonds ) {
    # 	    # First, checks if rotatable bond has fourth atom produce dihedral
    # 	    # angle. It is done by looking at atom connections - if rotatable
    # 	    # bond ends with terminal atom, then this bond is excluded.
    	    if( scalar( @{ $residue_site->
			       {$uniq_rotatable_bonds{$angle_name}->[1]}
			       {"connections"} } ) < 2 ){ next; }

    	    # Chooses proper atom ids for calculating dihedral angles.
    	    my $second_atom_id = $uniq_rotatable_bonds{$angle_name}->[0];
    	    my $third_atom_id = $uniq_rotatable_bonds{$angle_name}->[1];
    	    my @second_connections = # Second atom connections, except third.
    		grep { $_ ne $third_atom_id }
    	        @{ $residue_site->{$second_atom_id}{"connections"} };
    	    my $first_atom_name =
    		sort_by_priority(
    		filter( { "atom_site" => $residue_site,
    			  "include" => { "id" => \@second_connections },
    			  "data" => [ "label_atom_id" ],
    			  "is_list" => 1 } ) )->[0];
    	    my @third_connections = # Third atom connections, except second.
    		grep { $_ ne $second_atom_id }
    	        @{ $residue_site->{$third_atom_id}{"connections"} };
    	    my $fourth_atom_name =
    		sort_by_priority(
    		filter( { "atom_site" => $residue_site,
    			  "include" => { "id" => \@third_connections },
    			  "data" => [ "label_atom_id" ],
    			  "is_list" => 1 } ) )->[0];

    	    # Extracts coordinates for dihedral angle calculations.
    	    my $first_atom_coord =
    		filter(
    		    { "atom_site" => $residue_site,
    		      "include" => { "label_atom_id" => [ $first_atom_name ] },
    		      "data" => [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		      "is_list" => 1 } );
    	    my $second_atom_coord =
    		[ $residue_site->{$second_atom_id}{"Cartn_x"},
    		  $residue_site->{$second_atom_id}{"Cartn_y"},
    		  $residue_site->{$second_atom_id}{"Cartn_z"} ];
    	    my $third_atom_coord =
    		[ $residue_site->{$third_atom_id}{"Cartn_x"},
    		  $residue_site->{$third_atom_id}{"Cartn_y"},
    		  $residue_site->{$third_atom_id}{"Cartn_z"} ];
    	    my $fourth_atom_coord =
    		filter(
    		    { "atom_site" => $residue_site,
    		      "include" => { "label_atom_id" => [ $fourth_atom_name ] },
    		      "data" => [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    		      "is_list" => 1 } );

    	    $angle_values{$angle_name} =
    		dihedral_angle( [ $first_atom_coord,
    				  $second_atom_coord,
    				  $third_atom_coord,
    				  $fourth_atom_coord ] );
    	}

    	%{ $residue_angles{$residue_id} } = %angle_values;
    }

    return \%residue_angles;
}

#
# Calculates root-mean-square deviation of two same-length sets.
# Input:
#     ${first, second}_set - two equal by length sets of cartesian coordinates
#     of points.
# Output:
#     $rmsd - root-mean-square deviation.
#

sub rmsd
{
    my ( $first_set, $second_set ) = @_; # Order of set is not important.

    my $rmsd = 0;

    # Sums up sqaured differences of coordinates.
    for( my $i = 0; $i <= $#{ $first_set }; $i++ ) {
    	$rmsd += ( $first_set->[$i][0] - $second_set->[$i][0] )**2
    	       + ( $first_set->[$i][1] - $second_set->[$i][1] )**2
    	       + ( $first_set->[$i][2] - $second_set->[$i][2] )**2;
    }

    # Devides by the number of member of the set.
    $rmsd = $rmsd / scalar( @{ $first_set } );

    return sqrt( $rmsd );
}

1;
