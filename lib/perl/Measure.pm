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
use List::MoreUtils qw(uniq);

use lib "./";
use CifParser qw( filter_atoms
                  select_atom_data );
use LinearAlgebra qw( matrix_sub
                      vector_cross );
use LoadParams qw( rotatable_bonds );

my $parameter_file = "../../parameters/rotatable_bonds.csv";
use Data::Dumper;
# ----------------------------- Molecule parameters --------------------------- #

#
# Parameters.
#

my %ROTATABLE_BONDS = %{ rotatable_bonds( $parameter_file ) };

#
# Calculates various parameters that describe molecule or atoms, such as, bond
# length, dihedral angle, torsion angle, RMSD and etc.
#

#
# Calculates bond length of given two atoms.
# Input  (2 arg): matrices of x,y,z coordinates of two atoms.
# Output (1 arg): length of the bond in angstroms.
#

sub bond_length
{
    my @atom_coord = @_;

    my $bond_length =
	sqrt( ( $atom_coord[1][0] - $atom_coord[0][0] )**2
	    + ( $atom_coord[1][1] - $atom_coord[0][1] )**2
	    + ( $atom_coord[1][2] - $atom_coord[0][2] )**2 );

    return $bond_length;
}

#
# Calculates angle between three atoms.
# Input  (3 arg): matrices of x,y,z coordinates of three atoms.
# Output (1 arg): angle in radians.
#

sub bond_angle
{
    my $atom_coord = shift;
    my @atom_coord = @$atom_coord;

    my $bond_angle;

    # Angle between three atoms (in radians) in 3-D space can be calculated by
    # the formula:
    #                            ->   ->      ->         ->
    #            theta = arccos( AB * BC / || AB || * || BC || )

    # This formula is applied to atoms where vectors are the substraction of
    # coordinates of two atoms. Suppose, one of the side atom is A, B - middle
    # and C - remaining atom. Order of side atoms is irrelevant.
    my @vector_ab = ( $atom_coord[0][0] - $atom_coord[1][0],
		      $atom_coord[0][1] - $atom_coord[1][1],
		      $atom_coord[0][2] - $atom_coord[1][2] );
    my @vector_bc = ( $atom_coord[2][0] - $atom_coord[1][0],
		      $atom_coord[2][1] - $atom_coord[1][1],
		      $atom_coord[2][2] - $atom_coord[1][2] );

    my $vector_product = $vector_ab[0] * $vector_bc[0] +
	               + $vector_ab[1] * $vector_bc[1]
	               + $vector_ab[2] * $vector_bc[2];

    my $length_ab = sqrt( $vector_ab[0]**2
			+ $vector_ab[1]**2
			+ $vector_ab[2]**2 );
    my $length_bc = sqrt( $vector_bc[0]**2
			+ $vector_bc[1]**2
			+ $vector_bc[2]**2 );

    $bond_angle = acos( $vector_product / ( $length_ab * $length_bc ) );

    return $bond_angle;
}

#
# Calculates dihedral angle of four given atoms.
# Input  (4 arg): matrices of x,y,z coordinates of four atoms.
# Output (1 arg): dihedral angle in radians.
#

sub dihedral_angle
{
    my $atom_coord = shift;
    my @atom_coord = @$atom_coord;

    my $dihedral_angle;

    #                  -> ->    ->
    # Creates vectors: a, b and c, that are translated to global reference frame.
    # Picture of vectors:
    #                                   ->  O ->
    #                                   b  /  c
    #                              -> CA---C
    #                              a /
    #                               N
    my $vector_a = matrix_sub( [ $atom_coord[1] ], [ $atom_coord[0] ] );
    my $vector_b = matrix_sub( [ $atom_coord[2] ], [ $atom_coord[1] ] );
    my $vector_c = matrix_sub( [ $atom_coord[3] ], [ $atom_coord[2] ] );

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
    my @normal_vector_ab = map { $_ / $vector_length_ab } @$vector_cross_ab;
    my @normal_vector_bc = map { $_ / $vector_length_bc } @$vector_cross_bc;

    # Finishes orthonormal frame from normal vector ab, vector b and its cross
    # product.
    my $vector_length_b = sqrt( $vector_b->[0][0]**2
    			      + $vector_b->[0][1]**2
    			      + $vector_b->[0][2]**2 );
    my @normal_vector_b = map { $_ / $vector_length_b } @{ $vector_b->[0] };

    my @orthonormal_cross = vector_cross( \@normal_vector_ab, \@normal_vector_b );

    # Using orthonormal frame, projections from vector a and c are
    # generated and angle calculated.
    # TODO: check, if angle sign is properly assigned.
    $dihedral_angle =
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
# Input  (1 arg): atom site data structure.
# Output (1 arg): list of lists of dihedral angles.
#

sub all_dihedral
{
    my ( $atom_site ) = @_;


    # Selects specified residue(s) id(s).
    # TODO: code is similar to rotation_only in SidechainModels package. Possible
    # duplication.
    my $target_resi_id =
	&select_atom_data( [ "label_seq_id" ], $atom_site );
    my @target_resi_id = map { $_->[0] } @{ $target_resi_id };

    my @target_atom_id;

    my $residue_name;

    my $atom_type;
    my $mid_atom_type;
    my $up_atom_type;

    my $mid_atom_coord;
    my $up_atom_coord;
    my $side_atom_coord;

    my @rotatable_bonds;

    for my $resi_id ( uniq( @target_resi_id ) ) {

	print Dumper $ROTATABLE_BONDS{"SER"};

	    # my $angle_symbol; # Because side chain might have multiple
	    #                   # rotatable bonds, there must be distinct
	    #                   # symbols for different dihedral angles.

	    # for( my $i = 0; $i < scalar( @rotatable_bonds ); $i++ ) {
	    # 	$angle_symbol = "chi${i}";
	    # 	$mid_atom_type = $rotatable_bonds[$i][0];
	    # 	$mid_atom_type =~ s/\s//g;
	    # 	$up_atom_type = $rotatable_bonds[$i][1];
	    # 	$up_atom_type =~ s/\s//g;

    	#     # Information about side atom is stored in rotatable bonds array,
    	#     # except for CA atom.
    	#     if( $mid_atom_type eq "CA" ) {
    	# 	$side_atom_coord =
    	# 	    &select_atom_data(
    	# 	    [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    	# 	    &filter_atoms(
    	# 		{ "label_atom_id" => [ "N" ] },
    	# 		$atom_site ) );
    	#     } else {
    	# 	$side_atom_coord =
    	# 	    &select_atom_data(
    	# 	    [ "Cartn_x", "Cartn_y", "Cartn_z" ],
    	# 	    &filter_atoms(
    	# 		{ "label_atom_id" => [ $rotatable_bonds[$i-1][1] ] },
    	# 		$atom_site ) );
    	#     }

    	#     $mid_atom_coord =
    	# 	&select_atom_data(
    	# 	[ "Cartn_x", "Cartn_y", "Cartn_z" ],
    	# 	&filter_atoms(
    	# 	    { "label_atom_id" => [ $mid_atom_type ] },
    	# 	    $atom_site ) );

    	#     $up_atom_coord =
    	# 	&select_atom_data(
    	# 	[ "Cartn_x", "Cartn_y", "Cartn_z" ],
    	# 	&filter_atoms(
    	# 	    { "label_atom_id" => [ $up_atom_type ] },
    	# 	    $atom_site ) );

	    # }
    }
}


#
# Calculates root-mean-square deviation of two same-length sets.
# Input  (2 arg): two equal by length sets of cartesian coordinates of points.
# Output (1 arg): root-mean-square deviation.
#

sub rmsd
{
    my ( $first_set, $second_set ) = @_; # Order of set is not important.

    my $rmsd = 0; # TODO: check, if 0 will produce future bugs.

    # Sums up sqaured differences of coordinates.
    for( my $i = 0; $i <= $#{ $first_set }; $i++ ) {
    	$rmsd += ( $first_set->[$i][0] - $second_set->[$i][0] )**2
    	       + ( $first_set->[$i][1] - $second_set->[$i][1] )**2
    	       + ( $first_set->[$i][2] - $second_set->[$i][2] )**2;
    }

    # Devides by the number of member of the set.
    $rmsd = $rmsd / scalar( @$first_set );

    return sqrt( $rmsd ); # Takes square root.
}

1;
