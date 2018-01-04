package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( add_hydrogens
                     generate_library
                     generate_pseudo
                     generate_rotamer );

use List::Util qw( max );
use Math::Trig qw( acos );

use AtomInteractions qw( potential );
use AtomProperties qw( %ATOMS );
use Combinatorics qw( permutation );
use ConnectAtoms qw( connect_atoms );
use LinearAlgebra qw( find_euler_angles
                      flatten
                      mult_matrix_product
                      matrix_sum
                      matrix_sub
                      pi
                      reshape
                      scalar_multipl
                      switch_ref_frame
                      translation
                      vector_cross
                      vector_length
                      x_axis_rotation
                      y_axis_rotation
                      z_axis_rotation );
use Measure qw( all_dihedral
                dihedral_angle
                bond_angle
                bond_length );
use MoleculeProperties qw( %ROTATABLE_BONDS
                           %HYBRIDIZATION
                           %HYDROGEN_NAMES );
use PDBxParser qw( atom_data_with_id
                   create_pdbx_entry
                   filter_atoms
                   select_atom_data );
use Sampling qw( sample_angles );
use Data::Dumper;
# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates pseudo-atoms from side chain models that are written in equation
# form.
# Input:
#     $atom_site - atom site data structure (see PDBxParser). Must have any
#     sidechain model function applied to it (see SidechainModels.pm).
#     $atom_specifier - hash of hashes for selecting atoms by attributes (see
#     PDBxParser).
#     $angle_values - hash of arrays that describe possible values of dihedral
#     angles.
#     Ex.: { "chi0" => [ 0, 0.4, 1.5, 2.0 ],
#            "chi1" => [ 0, 2 ] }
# Output:
#     $pseudo_atom_site - atom site data structure for pseudo-atoms.
#

sub generate_pseudo
{
    my ( $atom_site, $atom_specifier, $angle_values ) = @_;

    my %pseudo_atom_site;

    my @atom_ids =
    	map { $_->[0]->[0] }
    	select_atom_data(
    	filter_atoms( $atom_site, $atom_specifier ),
    	[ "id" ] );
    my $last_atom_id = max( keys %{ $atom_site } );

    for my $atom_id ( @atom_ids ) {
    	# Calculates current dihedral angles of rotatable bonds. Will be used
    	# for reseting dihedral angles to 0 degree angle.
    	my $residue_id = $atom_site->{"$atom_id"}{"label_seq_id"};

    	my %angles =
	    %{ all_dihedral(
		   filter_atoms( $atom_site,
				 { "label_seq_id" => [ $residue_id ],
				   "label_alt_id" => [ "." ] } ) ) };

    	# Iterates through combinations of angles and evaluates conformational
    	# model.
    	my @angle_names = sort( { $a cmp $b } keys %{ $angle_values } );
    	my @angle_values;
    	for my $angle_name ( @angle_names ) {
    	    push( @angle_values,
    		  [ map
    		    { $_ - $angles{"$residue_id"}{"$angle_name"} }
    		    @{ $angle_values->{"$angle_name"} } ] );
    	}

    	my $conformation = $atom_site->{"$atom_id"}{"conformation"};
    	for my $angle_comb ( # Abreviation of angle combinations.
    	    @{ permutation( scalar( @angle_names ), [], \@angle_values, [] ) } ){
    	    my %angle_values =
    		map { ( $angle_names[$_] => $angle_comb->[$_] ) } 0..$#angle_names;
    	    # Evaluates matrices.
    	    my ( $transf_atom_coord ) =
    		@{ mult_matrix_product( $conformation, \%angle_values ) };

    	    # Adds necessary PDBx entries to pseudo atom site.
    	    $last_atom_id++;
    	    create_pdbx_entry(
    	    	{ "atom_site" => \%pseudo_atom_site,
    	    	  "id" => $last_atom_id,
    	    	  "type_symbol" => $atom_site->{$atom_id}{"type_symbol"},
    	    	  "label_atom_id" => $atom_site->{$atom_id}{"label_atom_id"},
    	    	  "label_alt_id" => "1",
    	    	  "label_comp_id" => $atom_site->{$atom_id}{"label_comp_id"},
    	    	  "label_asym_id" => $atom_site->{$atom_id}{"label_asym_id"},
    	    	  "label_entity_id" => $atom_site->{$atom_id}{"label_entity_id"},
    	    	  "label_seq_id" => $residue_id,
    	    	  "cartn_x" => sprintf( "%.3f", $transf_atom_coord->[0][0] ),
    	    	  "cartn_y" => sprintf( "%.3f", $transf_atom_coord->[1][0] ),
    	    	  "cartn_z" => sprintf( "%.3f", $transf_atom_coord->[2][0] ),
    		  "auth_seq_id" => $atom_site->{$atom_id}{"auth_seq_id"} } );
    	    # Adds information about used dihedral angles.
    	    $pseudo_atom_site{$last_atom_id}{"dihedral_angles"} =
    	    	\%angle_values;
    	    # Adds additional pseudo-atom flag for future filtering.
    	    $pseudo_atom_site{$last_atom_id}{"is_pseudo_atom"} = 1;
    	    # Adds atom id that pseudo atoms was made of.
    	    $pseudo_atom_site{$last_atom_id}{"origin_atom_id"} = $atom_id;
    	}
    }

    return \%pseudo_atom_site;
}

#
# Generates rotamers according to given angle values.
# Input:
#     $residue_atom_site - atom site data structure (see PDBxParser).
#     $angle_values - name and value of angles in hash form.
# Output:
#     %generated_rotamers - atom site data structure with additional
#     rotamer data.
#

sub generate_rotamer
{
    my ( $atom_site, $angle_values ) = @_;

    # Extracts residue name from residue id. Only one ID can be parsed.
    # TODO: maybe should consider generating rotamers for multiple residues.
    my @residue_ids = keys %{ $angle_values };
    my $residue_id = $residue_ids[0];
    my $residue_name =
    	select_atom_data(
	filter_atoms( $atom_site,
	{ "label_seq_id" => [ $residue_id ] } ),
	[ "label_comp_id" ] )->[0][0];

    my $atom_names =
    	select_atom_data(
    	filter_atoms( $atom_site,
    	{ "label_seq_id" => [ $residue_id ],
	  "label_alt_id" => [ "." ] } ),
	[ "label_atom_id" ] );
    my @atom_names = map { $_->[0] } @{ $atom_names };

    # Iterates through every atom of certain residue name and rotates to
    # specified dihedral angle.
    my %rotamer_atom_site;

    for my $atom_name ( @atom_names ) {
    	# Defines dihedral angles by %ROTATABLE_BONDS description and specified
    	# angles in $angles.
    	my %angles;
    	if( defined $ROTATABLE_BONDS{"$residue_name"}{"$atom_name"} ) {
    	    my $atom_id =
    	    	select_atom_data(
    		filter_atoms( $atom_site,
    		{ "label_atom_id" => [ $atom_name ],
		  "label_alt_id"  => [ "." ] } ),
    	    	[ "id" ] )->[0][0];

    	    for my $angle_id
    		( 0..scalar( @{ $ROTATABLE_BONDS{"$residue_name"}
    				                {"$atom_name"} } ) - 2 ) {
    		    $angles{"chi$angle_id"} =
    			[ $angle_values->{"$residue_id"}{"chi$angle_id"} ];
    	    }

    	    %rotamer_atom_site =
    		( %rotamer_atom_site,
    		  %{ generate_pseudo( { %{ $atom_site }, %rotamer_atom_site },
    				      { "id" => [ $atom_id ] },
    				      \%angles ) } );
    	}
    }

    return \%rotamer_atom_site;
}

#
# Generates rotamer libraries by specified arguments that include atom
# movements and interactions between atoms.
# Input:
#     $atom_site - atom site data structure (see PDBxParser).
#     $residue_ids - array of residue ids.
#     $small_angle - angle by which rotation is made.
#     $conformations - possible sidechain movements described by sidechain
#     modeling
#     functions in SidechainModels.pm.
#     $interactions - interaction models described by functions in
#     AtomInteractions.pm
# Output:
#     %library_atom_site - atom site data structure with additional data
#

sub generate_library
{
    my ( $args ) = @_;
    my $atom_site = $args->{"atom_site"};
    my $residue_ids = $args->{"residue_ids"};
    my $small_angle = $args->{"small_angle"};
    my $conformations = $args->{"conformations"};
    my $interactions = $args->{"interactions"};
    my $cutoff = $args->{"cutoff"};

    my %library_atom_site;

    # Generates comformational models before checking for clashes/interactions.
    $conformations->( $atom_site );

    my @sampled_angles =
	map { [ $_ ] } @{ sample_angles( [ [ 0, 2 * pi() ] ], $small_angle ) };

    for my $residue_id ( @{ $residue_ids } ) {
    	my $residue_site =
    	    filter_atoms( $atom_site, { "label_seq_id" => [ $residue_id ] } );
    	my $residue_name =
    	    select_atom_data( $residue_site, [ "label_comp_id" ] )->[0][0];
	my $atom_names =
	    select_atom_data( $residue_site, [ "label_atom_id" ] );

	# Creates a list of atoms that depend on the rotation of bonds and also
	# exist in current residue.
	my @sorted_names;
	for my $atom_name ( @{ $atom_names } ) {
	    push( @sorted_names, $atom_name->[0] )
		if exists $ROTATABLE_BONDS{"$residue_name"}{$atom_name->[0]};
	}

    	# Sorts atom names by the quantity of rotatable bonds described in
    	# ROTATABLE_BONDS.
    	@sorted_names =
    	    sort{ scalar( @{ $ROTATABLE_BONDS{"$residue_name"}{$a} } )
    	      cmp scalar( @{ $ROTATABLE_BONDS{"$residue_name"}{$b} } ) }
	    @sorted_names;

	# Ignores movable side chain atoms so, iteractions between pseudo atoms
	# and backbone could be analyzed properly.
	my %interaction_site = %{ $atom_site };
	my @removable_atom_ids =
	    map { $_->[0] }
	    @{ select_atom_data(
	       filter_atoms( $residue_site,
	       { "label_atom_id" => \@sorted_names } ),
	       [ "id" ] ) };
	for my $atom_id ( @removable_atom_ids ) {
	    delete $interaction_site{"$atom_id"}
	}

    	# Iterates through sorted atoms and tries to detect interactions.
	my @allowed_angles; # comb - combinations.

    	# my %current_angles;
    	for my $atom_name ( @sorted_names ) {
    	    my $current_atom_site =
    	    	filter_atoms( $residue_site,
    	    		      { "label_atom_id" => [ "$atom_name" ] } );
    	    my $current_atom_id =
    	    	select_atom_data( $current_atom_site, [ "id" ] )->[0][0];

	    my $angle_count =
		scalar( @{ $ROTATABLE_BONDS{"$residue_name"}
			                   {"$atom_name"} } ) - 2;

	    # TODO: look for cases, when all atoms produce clashes.
	    my @current_angles;
	    if( ! @allowed_angles ) {
		# If no angles are present, allows all angles.
		@current_angles = @sampled_angles;
	    } else {
		my $allowed_angle_count = scalar( $allowed_angles[0] );
		if( $angle_count == $allowed_angle_count ) {
		    @current_angles = @allowed_angles;
		} else {
		    @current_angles =
		    	@{ permutation( 2, [], [ \@allowed_angles,
		    				 \@sampled_angles ], [ ] ) };
		    # Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
		    @current_angles =
		    	map { [ @{ $_->[0] }, @{ $_->[1] } ] } @current_angles;
		}
	    }

	    undef @allowed_angles; # Undefying to store increasing dihedral
	                           # angles.

	    # Iterates through allowed angles, checks clashes and returns
	    # unhindered angles to @allowed_angle_comb.
	    for my $angles ( @current_angles ) {
		my %angles =
		    map { ( "chi$_" => [ $angles->[$_] ] ) } 0..$angle_count;

		my $pseudo_atom_site =
		    generate_pseudo( $atom_site,
				     { "id" => [ "$current_atom_id" ] },
				     \%angles );
		my $pseudo_atom_id =
		    select_atom_data( $pseudo_atom_site, [ "id" ] )->[0][0];

		potential( { %interaction_site,
			     %{ $pseudo_atom_site } },
			   $interactions,
			   $cutoff,
			   { "id" => [ $pseudo_atom_id ] },
			   { "label_atom_id" => [ "N", "CA", "C", "O" ] } );

		if( $pseudo_atom_site->{"$pseudo_atom_id"}
		                       {"potential_energy"} <= $cutoff ) {
		    push( @allowed_angles, $angles );
		}
	    }

	    die "No possible rotamer solutions were detected."
	    	if scalar( @allowed_angles ) == 0;
	}

	# Generates final rotamers.
	for my $angles ( @allowed_angles ) {
	    my %angles =
		( "$residue_id" => { map { ( "chi$_" => $angles->[$_] ) }
				     ( 0..$#{ $angles } ) } );
	    # TODO: do not forget to check clashes among atoms inside rotamer.
	    %library_atom_site =
	    	( %library_atom_site,
	    	  %{ generate_rotamer( { %library_atom_site,
	    				 %{ $atom_site } },
				       \%angles ) } );
	}
    }

    return \%library_atom_site;
}

sub add_hydrogens
{
    my ( $atom_site ) = @_;

    my %atom_site = %{ $atom_site };
    connect_atoms( \%atom_site );

    my %hydrogen_site;
    my $last_atom_id = max( keys %{ $atom_site } );

    for my $atom_id ( sort { $a <=> $b } keys %atom_site ) {
    	my $atom_type = $atom_site{$atom_id}{"type_symbol"};
    	my $atom_name = $atom_site{$atom_id}{"label_atom_id"};

	my $residue_name = $atom_site{$atom_id}{"label_comp_id"};
	my $residue_id = $atom_site{$atom_id}{"label_seq_id"};
	my $residue_site =
	    filter_atoms( $atom_site, { "label_seq_id" => [ $residue_id ] } );
	my %residue_coord =
	    %{ atom_data_with_id( $residue_site,
				  [ "Cartn_x", "Cartn_y", "Cartn_z" ] ) };

    	my $hydrogen_names = $HYDROGEN_NAMES{$residue_name}{$atom_name};

    	if( ! $hydrogen_names ) { next; }; # Exits early if there should be no
	                                   # hydrogens connected to the atom.

    	my @connection_ids = @{ $atom_site{"$atom_id"}{"connections"} };
    	my @connection_names =
    	    map { $atom_site{"$_"}{"label_atom_id"} } @connection_ids;

    	my $hybridization = $HYBRIDIZATION{$residue_name}{$atom_name};
    	my $lone_pair_count = $ATOMS{$atom_type}{"lone_pairs"};

    	# Decides how many and what hydrogens should be added according to the
    	# quantity of bonds and hydrogen atoms that should be connected to the
	# target atom.
    	my @missing_hydrogens;
    	for my $hydrogen_name ( @{ $hydrogen_names } ) {
    	    if( ! grep { /$hydrogen_name/ } @connection_names ) {
    	    	push( @missing_hydrogens, $hydrogen_name );
    	    }
    	}

    	# Exits early if there are no spare hydrogens to add.
    	if( scalar( @missing_hydrogens ) == 0 ) { next; };

    	#            sp3                       sp2               sp
    	#
    	#            Up(2)                     Up(2)	         Up(2)
    	# z          |		               |	         |
    	# |_y      Middle(1) __ Right(3)     Middle(1)	       Middle(1)
    	# /         / \		              / \	         |
    	# x    Left(4) Back(5)	         Left(4) Right(3)        Down(3)
    	#
    	# Depending on hybridization and present bond connections, adds missing
    	# hydrogens.
    	# TODO: in the future, should adjust sp3, sp2 angles according to
    	# experimental data, not model.
    	my %hydrogen_coord = map { $_ => undef } @missing_hydrogens;
    	# my $transf_matrix;

    	if( $hybridization eq "sp3" ) {
	    # Because bond length depends on hybridization of atoms, bond length
	    # is calculated for each differently hybridized atoms.
    	    my $bond_length =
    		$ATOMS{$atom_type}{"covalent_radius"}{"length"}[0]
    	      + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

    	    if( scalar( @connection_ids ) == 3 ) {
		my ( $up_atom_coord,
		     $mid_atom_coord,
		     $left_atom_coord,
		     $right_atom_coord ) =
			 ( $residue_coord{$connection_ids[0]},
			   $residue_coord{$atom_id},
			   $residue_coord{$connection_ids[1]},
			   $residue_coord{$connection_ids[2]} );

    	    	# Theoretically optimal angles should be so, the distance of
    	    	# atoms would be furthest from each other. This strategy could
    	    	# be achieved by imagining each bond as vector of the length 1.
    	    	# Then the sum of all vectors would be 0, if angles are equal.
    	    	#
    	    	#                        ->  ->  ->  ->
    	    	#                        A + B + C + D = 0
    	    	#
    	    	# In that case, the square of sum also should be equal 0.
    	    	#
    	    	#                        ->  ->  ->  ->
    	    	#                      ( A + B + C + D ) ^ 2 = 0
    	    	#
    	    	#       ->  ->      ->  ->      ->  ->      ->  ->      ->  ->
    	    	#   2 * A * B + 2 * A * C + 2 * A * D + 2 * B * C + 2 * B * D +
    	    	#       ->  ->   ->      ->      ->      ->
    	    	# + 2 * C + D  + A ^ 2 + B ^ 2 + C ^ 2 + D ^ 2 = 0
    	    	#
    	    	# Because length of each vector is equal to 1, then dot product
    	    	# is equal to cos(alpha), where all angles between bonds are
    	    	# equal. If there were no given angles, alpha should be 109.5
    	    	# degrees.
    	    	#
    	    	#                   alpha = arccos( - 1 / 3 )
    	    	#
    	    	# However, there is a restriction of given angle. And the
    	    	# calculation changes:
    	    	#
    	    	#        alpha = arccos( ( - 4 - 2 * cos( beta )
		#                              - 2 * cos( gamma )
		#                              - 2 * cos( delta ) ) / 6 )
    	    	#
    	    	# where beta is the given angle.

		# Calculates all bond angles between atoms connected to the
		# middle atom.
		my $up_mid_right_angle =
		    bond_angle( [ $up_atom_coord,
				  $mid_atom_coord,
				  $right_atom_coord ] );
		my $up_mid_left_angle =
		    bond_angle( [ $up_atom_coord,
				  $mid_atom_coord,
				  $left_atom_coord ] );
		my $right_mid_left_angle =
		    bond_angle( [ $up_atom_coord,
				  $mid_atom_coord,
				  $left_atom_coord ] );

		# Calculates what common angle between hydrogen and rest of the
		# atoms should be.
		my $hydrogen_angle =
		    acos( ( - 4
			    - 2 * cos( $up_mid_right_angle )
			    - 2 * cos( $up_mid_left_angle )
			    - 2 * cos( $right_mid_left_angle ) )
			  / 6 );

		# Determines dihedral angle between left and right atoms. Then
		# splits rest of the 2 * pi angle into two equal parts.
		my $dihedral_angle =
		    dihedral_angle( [ $left_atom_coord,
				      $up_atom_coord,
				      $mid_atom_coord,
				      $right_atom_coord ] );
		if( abs( $dihedral_angle ) < ( 3 * pi() / 4 ) ) {
		    if( $dihedral_angle < 0 ) {
			$dihedral_angle = ( 2 * pi() + $dihedral_angle ) / 2;
		    } else {
			$dihedral_angle = - ( 2 * pi() - $dihedral_angle ) / 2;
		    }
		} else {
		    if( $dihedral_angle < 0 ) {
			$dihedral_angle = $dihedral_angle / 2;
		    } else {
			$dihedral_angle = - $dihedral_angle / 2;
		    }
		}

		# Places hydrogen according to previously calculated angles.
    	    	my ( $transf_matrix ) =
    	    	    @{ switch_ref_frame(
			   $mid_atom_coord,
			   $up_atom_coord,
			   $left_atom_coord,
			   "global" ) };

    	    	( $hydrogen_coord{$missing_hydrogens[0]} ) =
    	    	    @{ mult_matrix_product(
    	    		   [ $transf_matrix,
    	    		     [ [ $bond_length
    	    			* cos( pi() / 2 - $dihedral_angle )
    	    		        * sin( $hydrogen_angle ) ],
    	    		       [ $bond_length
    	    		        * sin( pi() / 2 - $dihedral_angle )
    	    		        * sin( $hydrogen_angle ) ],
    	    		       [ $bond_length
    	    		        * cos( $hydrogen_angle ) ],
    	    		       [ 1 ] ] ] ) };
    	    } elsif( scalar( @connection_ids ) == 2 ) {
    	    	# Calculates current angle between atoms that are connected to
    	    	# target atom.
		my ( $up_atom_coord,
		     $mid_atom_coord,
		     $left_atom_coord ) =
			 ( $residue_coord{$connection_ids[0]},
			   $residue_coord{$atom_id},
			   $residue_coord{$connection_ids[1]} );

    	    	my $bond_angle =
    	    	    bond_angle( [ $up_atom_coord,
    	    			  $mid_atom_coord,
    	    			  $left_atom_coord ] );

    	    	# This time is only one defined angle. And the calculation
		# changes:
    	    	#
		#       acos( ( - 4 - 2 * cos( $bond_angle ) ) / 10 )
		#
    	    	my $hydrogen_angle =
    	    	    acos( ( - 4 - 2 * cos( $bond_angle ) ) / 10 );

    	    	# Generates transformation matrix for transfering atoms to local
    	    	# reference frame.
    	    	my ( $transf_matrix ) =
    	    	    @{ switch_ref_frame( $mid_atom_coord,
    	    				 $up_atom_coord,
    	    				 $left_atom_coord,
    	    				 'global' ) };

    	    	# Adds hydrogen first to both atoms that have 0 or 1 electron
    	    	# pairs.
    	    	if( scalar( @missing_hydrogens ) >= 1 ) {
    	    	    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    	    		@{ mult_matrix_product(
    	    		       [ $transf_matrix,
    	    			 [ [ $bond_length
    	    			   * cos( 7 * pi() / 6 )
    	    			   * sin( $hydrogen_angle ) ],
    	    			   [ $bond_length
    	    		           * sin( 7 * pi() / 6 )
    	    		           * sin( $hydrogen_angle )],
    	    			   [ $bond_length
    	    		           * cos( $hydrogen_angle )],
    	    			   [ 1 ] ] ] ) };
    	    	    shift @missing_hydrogens;
    	    	}

    	    	# Additional hydrogen is added only to the atom that has no
    	    	# electron pairs.
    	    	if( scalar( @missing_hydrogens ) == 1 ) {
    	    	    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    	    		@{ mult_matrix_product(
    	    		       [ $transf_matrix,
    	    			 [ [ $bond_length
    	    			   * cos( - 1 * pi() / 6 )
    	    			   * sin( $hydrogen_angle ) ],
    	    			   [ $bond_length
    	    		           * sin( - 1 * pi() / 6 )
    	    		           * sin( $hydrogen_angle ) ],
    	    			   [ $bond_length
    	    		           * cos( $hydrogen_angle ) ],
    	    			   [ 1 ] ] ] ) };
    	    	}

    	    } elsif( scalar( @connection_ids ) == 1 ) {
    	    	# Calculates current angle between atoms that are connected to
    	    	# target atom.
		my ( $up_atom_coord,
		     $mid_atom_coord,
		     $side_coord ) = # Coordinate that only will be used for
			             # defining a local reference frame.
			 ( $residue_coord{$connection_ids[0]},
			   $residue_coord{$atom_id},
			   [ $atom_site->{$atom_id}{"Cartn_x"},
			     $atom_site->{$atom_id}{"Cartn_y"} + 1,
			     $atom_site->{$atom_id}{"Cartn_z"} ] );

    	    	# Generates transformation matrix for transfering atoms to local
    	    	# reference frame.
    	    	my ( $transf_matrix ) =
    	    	    @{ switch_ref_frame( $mid_atom_coord,
    	    				 $up_atom_coord,
    	    				 $side_coord,
    	    				 'global' ) };

    	    	# Decreases bond angle, if lone pairs are present.
		# TODO: check if angle reduction is relevant in amino acid
		# structures.
    	    	my $bond_angle;
    	    	if( $lone_pair_count > 0 ) {
    	    	    $bond_angle =
    	    		( 109.5 - $lone_pair_count * 2.5 ) * pi() / 180;
    	    	} else {
    	    	    $bond_angle = 109.5 * pi() / 180;
    	    	}

    	    	# Adds hydrogens according to the quantity of lone pairs.
    	    	if( scalar( @missing_hydrogens ) >= 3 ) {
    	    	    ( $hydrogen_coord{$missing_hydrogens[0]} )=
    	    		@{ mult_matrix_product(
    	    		       [ $transf_matrix,
    	    			 [ [ $bond_length * sin( $bond_angle ) ],
    	    			   [ 0 ],
    	    			   [ $bond_length * cos( $bond_angle ) ],
    	    			   [ 1 ] ] ] ) };
    	    	    shift @missing_hydrogens;
    	    	}

    	    	if( scalar( @missing_hydrogens ) >= 2 ) {
    	    	    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    	    		@{ mult_matrix_product(
    	    		       [ $transf_matrix,
    	    			 [ [ $bond_length
    	    			   * cos( 2 * pi() / 3 )
    	    		           * sin( $bond_angle ) ],
    	    			   [ $bond_length
    	    		           * sin( 2 * pi() / 3 )
    	    		           * sin( $bond_angle ) ],
    	    			   [ $bond_length
    	    		           * cos( $bond_angle ) ],
    	    			   [ 1 ] ] ] ) };
    	    	    shift @missing_hydrogens;
    	    	}

    	    	if( scalar( @missing_hydrogens ) == 1 ) {
    	    	    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    	    		@{ mult_matrix_product(
    	    		       [ $transf_matrix,
    	    			 [ [ $bond_length
    	    			   * cos( 4 * pi() / 3 )
    	    			   * sin( $bond_angle ) ],
    	    			   [ $bond_length
    	    		           * sin( 4 * pi() / 3 )
    	    		           * sin( $bond_angle ) ],
    	    			   [ $bond_length
    	    		           * cos( $bond_angle ) ],
    	    			   [ 1 ] ] ] ) };
    	    	}
    	    }

    	} elsif( $hybridization eq "sp2" ) {
    	    my $bond_length =
    		$ATOMS{$atom_type}{"covalent_radius"}{"length"}[1]
    	      + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

    	    # Depending on quantity of atoms connections, adds hydrogens.
    	    if( scalar( @connection_ids ) == 2 ) {
    	    	# Calculates current angle between atoms that are connected to
    	    	# target atom.
		my ( $up_atom_coord,
		     $mid_atom_coord,
		     $left_atom_coord ) =
			 ( $residue_coord{$connection_ids[0]},
			   $residue_coord{$atom_id},
			   $residue_coord{$connection_ids[1]} );

    	    	# Generates transformation matrix for transfering atoms to local
    	    	# reference frame.
    		my ( $transf_matrix ) =
    		    @{ switch_ref_frame( $mid_atom_coord,
    					 $up_atom_coord,
    					 $left_atom_coord,
    					 'global' ) };

    		# Calculates angle between two bonds where target atom takes
    		# part. Then desirable hydrogen angle should be calculated by
    		# formula:
    		#
    		# angle = ( 360 - beta ) / 2
    		#
    		# where beta is angle between two given bonds.
    		my $bond_angle =
    		    ( 2 * pi() - bond_angle( [ $up_atom_coord,
    					       $mid_atom_coord,
    					       $left_atom_coord ] ) )
    		    / 2;

    		# Hydrogen is placed by placing hydrogen colinearly and then
    		# rotating according to bond angle.
    		( $hydrogen_coord{$missing_hydrogens[0]} ) =
    		    @{ mult_matrix_product(
    			   [ $transf_matrix,
    			     [ [ $bond_length
    			       * cos( - 0.5 * pi() )
    			       * sin( $bond_angle ) ],
    			       [ $bond_length
    			       * sin( - 0.5 * pi() )
    			       * sin( $bond_angle ) ],
    			       [ $bond_length
    			       * cos( $bond_angle ) ],
    			       [ 1 ] ] ] ) };

    	    } elsif( scalar( @connection_ids ) == 1 ) {
		my ( $up_atom_coord,
		     $mid_atom_coord,
		     $side_coord ) =
			 ( $residue_coord{$connection_ids[0]},
			   $residue_coord{$atom_id},
			   [ $atom_site->{$atom_id}{"Cartn_x"},
			     $atom_site->{$atom_id}{"Cartn_y"} + 1,
			     $atom_site->{$atom_id}{"Cartn_z"} ] );

    		# If terminal atom belongs to conjugated system, hydrogens are
    		# added not to violate rule where atoms should be in one plain.
    		my @second_neighbours =
    		    grep { ! /$atom_id/ }
    		    map { @{ $atom_site->{$_}{"connections"} } }
    		    @connection_ids;

    		for my $second_neighbour ( @second_neighbours ) {
    		    my $second_hybridization =
    			$HYBRIDIZATION{$residue_name}
    		                      {$atom_site->{$second_neighbour}
    				                   {"label_atom_id"}};
    		    if( $second_hybridization eq "sp2"
    		     || $second_hybridization eq "sp" ) {
    			$side_coord =
    			    [ $atom_site->{$second_neighbour}{"Cartn_x"},
    			      $atom_site->{$second_neighbour}{"Cartn_y"},
    			      $atom_site->{$second_neighbour}{"Cartn_z"} ];
    			last;
    		    }
    		}

    	    	# Generates transformation matrix for transfering atoms to local
    	    	# reference frame.
    		my ( $transf_matrix ) =
    		    @{ switch_ref_frame( $mid_atom_coord,
    					 $up_atom_coord,
    					 $side_coord,
    					 'global' ) };

    		if( scalar( @missing_hydrogens ) >= 1 ) {
    		    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    			@{ mult_matrix_product(
    			       [ $transf_matrix,
    				 [ [ $bond_length
    				   * cos( -0.5 * pi() )
    				   * sin( 120 * pi() / 180 ) ],
    				   [ $bond_length
    			           * sin( -0.5 * pi() )
    			           * sin( 120 * pi() / 180 ) ],
    				   [ $bond_length
    			           * cos( 120 * pi() / 180 ) ],
    				   [ 1 ] ] ] ) };
    		    shift @missing_hydrogens;
    		}

    		if( scalar( @missing_hydrogens ) == 1 ) {
    		    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    			@{ mult_matrix_product(
    			       [ $transf_matrix,
    				 [ [ $bond_length
    				   * cos( 0.5 * pi() )
    				   * sin( 120 * pi() / 180 ) ],
    				   [ $bond_length
    			           * sin( 0.5 * pi() )
    			           * sin( 120 * pi() / 180 ) ],
    				   [ $bond_length
    			           * cos( 120 * pi() / 180 ) ],
    				   [ 1 ] ] ] ) };
    		}
    	    }

    	} elsif( $hybridization eq "sp" ) {
    	    # Calculates length of bonds.
    	    my $bond_length =
    	    	$ATOMS{$atom_name}{"covalent_radius"}{"length"}[2]
    	      + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

	    my ( $up_atom_coord,
		 $mid_atom_coord,
		 $side_coord ) =
		     ( $residue_coord{$connection_ids[0]},
		       $residue_coord{$atom_id},
		       [ $atom_site->{$atom_id}{"Cartn_x"},
			 $atom_site->{$atom_id}{"Cartn_y"} + 1,
			 $atom_site->{$atom_id}{"Cartn_z"} ] );

	    # Generates transformation matrix for transfering atoms to local
	    # reference frame.
    	    my ( $transf_matrix ) =
    		@{ switch_ref_frame( $mid_atom_coord,
    				     $up_atom_coord,
    				     $side_coord,
    				     'global' ) };

    	    ( $hydrogen_coord{$missing_hydrogens[0]} ) =
    		@{ matrix_product(
    		       [ $transf_matrix,
    			 [ [ 0 ],
    			   [ 0 ],
    			   [ - $bond_length ],
    			   [ 1 ] ] ] ) };
    	}

    	# Each coordinate of atoms is transformed by transformation
    	# matrix and added to %hydrogen_site.
    	for my $hydrogen_name ( sort { $a cmp $b } keys %hydrogen_coord ) {
	    if( $hydrogen_coord{$hydrogen_name} ) {
	    # Adds necessary PDBx entries to pseudo atom site.
    	    $last_atom_id++;
	    create_pdbx_entry(
	    	{ "atom_site" => \%hydrogen_site,
	    	  "id" => $last_atom_id,
	    	  "type_symbol" => "H",
	    	  "label_atom_id" => $hydrogen_name,
	    	  "label_alt_id" => "1",
	    	  "label_comp_id" => $atom_site->{$atom_id}{"label_comp_id"},
	    	  "label_asym_id" => $atom_site->{$atom_id}{"label_asym_id"},
	    	  "label_entity_id" => $atom_site->{$atom_id}{"label_entity_id"},
	    	  "label_seq_id" => $residue_id,
	    	  "cartn_x" =>
		      sprintf( "%.3f", $hydrogen_coord{$hydrogen_name}->[0][0] ),
	    	  "cartn_y" =>
		      sprintf( "%.3f", $hydrogen_coord{$hydrogen_name}->[1][0] ),
	    	  "cartn_z" =>
		      sprintf( "%.3f", $hydrogen_coord{$hydrogen_name}->[2][0] ),
		  "auth_seq_id" => $atom_site->{$atom_id}{"auth_seq_id"}
		} );
    	    # Adds additional pseudo-atom flag for future filtering.
    	    $hydrogen_site{$last_atom_id}{"is_pseudo_atom"} = 1;
	    # Adds atom id that pseudo atoms was made of.
    	    $hydrogen_site{$last_atom_id}{"origin_atom_id"} = $atom_id;
	    }
    	}
    }

    return \%hydrogen_site;
}

1;
