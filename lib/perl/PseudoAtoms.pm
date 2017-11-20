package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( add_hydrogens
                     generate_library
                     generate_pseudo
                     generate_rotamer );

use List::Util qw( max );

use lib qw( ./ );
use AtomInteractions qw( potential );
use AtomProperties qw( %ATOMS );
use Combinatorics qw( permutation );
use ConnectAtoms qw( connect_atoms );
use LinearAlgebra qw( evaluate_matrix
                      find_euler_angles
                      matrix_product
                      matrix_sum
                      pi
                      scalar_multipl
                      switch_ref_frame
                      translation
                      vectorize
                      x_axis_rotation
                      y_axis_rotation
                      z_axis_rotation );
use Measure qw( all_dihedral
                bond_length );
use MoleculeProperties qw( %ROTATABLE_BONDS
                           %HYBRIDIZATION );
use PDBxParser qw( create_pdbx_entry
                   filter_atoms
                   select_atom_data );
use Sampling qw( sample_angles );

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
			       "label_alt_id" => [ "." ]} ) ) };

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
    		map { $angle_names[$_], $angle_comb->[$_] } 0..$#angle_names;
    	    # Converts matrices to GiNaC compatable format and evaluates them.
    	    my $transf_atom_coord =
		evaluate_matrix( $conformation, \%angle_values );

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
		  "cartn_z" => sprintf( "%.3f", $transf_atom_coord->[2][0] ) } );
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
#     $conformations - possible sidechain movements described by sidechain modeling
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

    	# Sorts atom names by the quantity of rotatable bonds described in
    	# ROTATABLE_BONDS.
    	my @sorted_names =
    	    sort{ scalar( @{ $ROTATABLE_BONDS{"$residue_name"}{$a} } )
    	      cmp scalar( @{ $ROTATABLE_BONDS{"$residue_name"}{$b} } ) }
	    keys %{ $ROTATABLE_BONDS{"$residue_name"} };

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
	    if( not @allowed_angles ) {
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

    # TODO: might not work with single atoms. Look into it. Also, would be faster
    # if heavy atoms without any possible addition of hydrogens are ignored.
    my %hydrogen_site;
    my $last_atom_id = max( keys %{ $atom_site } );

    for my $atom_id ( keys %atom_site ) {
	my $residue_name = $atom_site{$atom_id}{"label_comp_id"};
	my $atom_type = $atom_site{$atom_id}{"type_symbol"};
	my $atom_name = $atom_site{$atom_id}{"label_atom_id"};

	my @connection_ids = @{ $atom_site{"$atom_id"}{"connections"} };
	my @connection_atoms =
	    map { $atom_site{"$_"}{"type_symbol"} } @connection_ids;

	my $hybridization = $HYBRIDIZATION{$residue_name}{$atom_name};
	my $lone_pair_count = $ATOMS{$atom_type}{"lone_pairs"};
	my @lone_pair_types =
	    map { "." } 1..$lone_pair_count if $lone_pair_count > 0;

	my $hydrogen_count = 4 - scalar( @connection_ids ) - $lone_pair_count;
	my @hydrogen_types =
	    map { "H" } 1..$hydrogen_count if $hydrogen_count > 0;
	my @hydrogen_names;
	for my $hydrogen_id ( 1..$hydrogen_count ) {
	    ( my $hydrogen_name = $atom_name ) =~
		s/$atom_type(.?)/H$1$hydrogen_id/;
	    push( @hydrogen_names, $hydrogen_name );
	}

	if( $hybridization eq "sp3" ) {
	    my $tetrahedron_coord =
		sp3_tetrahedron( $atom_type,
				 @connection_atoms,
				 @hydrogen_types,
				 @lone_pair_types );

	    # Deals with sp3 hybrized atoms that has only one connection with
	    # another heavy atom.
	    if( scalar( @connection_ids ) == 1 ) {
		my $transf_matrix =
		    switch_ref_frame(
			[ $atom_site->{$atom_id}{"Cartn_x"},
			  $atom_site->{$atom_id}{"Cartn_y"},
			  $atom_site->{$atom_id}{"Cartn_z"} ],
			[ $atom_site->{$connection_ids[0]}{"Cartn_x"},
			  $atom_site->{$connection_ids[0]}{"Cartn_y"},
			  $atom_site->{$connection_ids[0]}{"Cartn_z"} ],
			[ $atom_site->{$atom_id}{"Cartn_x"} + 1,   # Arbitrary
			  $atom_site->{$atom_id}{"Cartn_y"} + 1,   # point.
			  $atom_site->{$atom_id}{"Cartn_z"} + 1 ], #
			"global" );

		# Adds hydrogens.
		for my $hydrogen_idx ( 0..$hydrogen_count-1 ) {
		    my $transf_atom_coord =
			matrix_product(
			    $transf_matrix,
			    vectorize( $tetrahedron_coord->[2+$hydrogen_idx] ) );

		    # Adds necessary PDBx entries to pseudo atom.
		    $last_atom_id++;
		    create_pdbx_entry(
			{ "atom_site" => \%hydrogen_site,
			  "id" => $last_atom_id,
			  "type_symbol" => "H",
			  "label_atom_id" => $hydrogen_names[$hydrogen_idx],
			  "label_alt_id" => ".",
			  "label_comp_id" =>
			      $atom_site->{$atom_id}{"label_comp_id"},
			  "label_asym_id" =>
			      $atom_site->{$atom_id}{"label_asym_id"},
			  "label_entity_id" =>
			      $atom_site->{$atom_id}{"label_entity_id"},
			  "label_seq_id" =>
			      $atom_site->{$atom_id}{"label_seq_id"},
			  "cartn_x" =>
			      sprintf( "%.3f", $transf_atom_coord->[0][0] ),
			  "cartn_y" =>
			      sprintf( "%.3f", $transf_atom_coord->[1][0] ),
			  "cartn_z" =>
			      sprintf( "%.3f", $transf_atom_coord->[2][0] ) } );
		    # Adds additional pseudo-atom flag for future filtering.
		    $hydrogen_site{$last_atom_id}{"is_pseudo_atom"} = 1;
		    # Adds atom id that pseudo atoms was made of.
		    $hydrogen_site{$last_atom_id}{"origin_atom_id"} = $atom_id;
		}
	    }
	}
	# } elsif( $hybridization eq "sp2" ) {
	#     my $hydrogen_count =
	# 	3 - scalar( @connection_ids ) - $ATOMS{$atom_type}{"lone_pairs"};
	#     for my $hydrogen_id ( 1..$hydrogen_count ) {
	# 	( my $hydrogen_name = $atom_name ) =~
	# 	    s/$atom_type(.?)/H$1$hydrogen_id/;
	#     }
	# } elsif( $hybridization eq "sp" ) {
	#     my $hydrogen_count =
	# 	2 - scalar( @connection_ids ) - $ATOMS{$atom_type}{"lone_pairs"};
	#     for my $hydrogen_id ( 1..$hydrogen_count ) {
	# 	( my $hydrogen_name = $atom_name ) =~
	# 	    s/$atom_type(.?)/H$1$hydrogen_id/;
	#     }
	# }
    }

    return \%hydrogen_site;
}

#                                       Up(2)
# z                                     |
# |_y                                   C(1) __ Right(3)
# /                                    / \
# x                               Left(4) Back(5)

sub sp3_tetrahedron
{
    my ( @atom_names ) = @_;

    # Places first atom in the origin of the global frame of reference.
    my @atom_coord = ( [ 0, 0, 0 ] );

    # Places second atom on z-axis by the length of the bond.
    my $bond_length_c_up =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[1]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord, [ 0, 0, $bond_length_c_up ] );

    # Places third atom in xz-plane that has 109.5 deg angle between 2-1-3 atoms.
    # TODO: in the future, should adjust sp3 angles according to experimental
    # angles, not model.
    my $bond_length_c_right =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[2]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord,
	  [ $bond_length_c_right * sin( 109.5 * pi() / 180 ),
	    0,
	    $bond_length_c_right * cos( 109.5 * pi() / 180 ) ] );

    # Places fourth atom in a position of third atom and rotates 109.5 deg.
    my $bond_length_c_left =
    	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[3]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord,
    	  [   $bond_length_c_left
	    * cos( 120 * pi() / 180 )
	    * sin( 109.5 * pi() / 180 ),
    	      $bond_length_c_left
	    * sin( 120 * pi() / 180 ) * sin( 109.5 * pi() / 180 ),
    	      $bond_length_c_left
	    * cos( 109.5 * pi() / 180 ) ] );

    # Places fifth atom in a position of third atom and rotates 109.5 deg.
    my $bond_length_c_back =
    	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[3]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord,
    	  [   $bond_length_c_back
	    * cos( 240 * pi() / 180 )
	    * sin( 109.5 * pi() / 180 ),
              $bond_length_c_back
	    * sin( 240 * pi() / 180 )
	    * sin( 109.5 * pi() / 180 ),
              $bond_length_c_back
	    * cos( 109.5 * pi() / 180 ) ] );

    return \@atom_coord;
}

#                                       Up(2)
# z                                     |
# |_y                                   C(1)
# /                                    / \
# x                               Left(4) Right(3)

sub sp2_triangle
{
    my ( @atom_names ) = @_;

    # Places first atom in the origin of the global frame of reference.
    my @atom_coord = ( [ 0, 0, 0 ] );

    # Places second atom on z-axis by the length of the bond.
    my $bond_length_c_up =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[1]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord, [ 0, 0, $bond_length_c_up ] );

    # TODO: in the future, should adjust sp2 angles according to experimental
    # angles, not model.
    # Places third atom in xz-plane that has 120 deg angle between 2-1-3 atoms.
    my $bond_length_c_right =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[2]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord,
	  [ $bond_length_c_right * sin( 120 * pi() / 180 ),
	    0,
	    $bond_length_c_right * cos( 120 * pi() / 180 ) ] );

    # Places fourth atom in xz-plane that has 120 deg angle between 2-1-4 atoms.
    my $bond_length_c_left =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[2]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord,
	  [ $bond_length_c_left * sin( 240 * pi() / 180 ),
	    0,
	    $bond_length_c_left * cos( 240 * pi() / 180 ) ] );

    return \@atom_coord;

}

#                                       Up(2)
# z                                     |
# |_y                                   C(1)
# /                                     |
# x                                     Down(3)

sub sp_line
{
    my ( @atom_names ) = @_;

    # Places first atom in the origin of the global frame of reference.
    my @atom_coord = ( [ 0, 0, 0 ] );

    # Places second atom on z-axis by the length of the bond.
    my $bond_length_c_up =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[1]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord, [ 0, 0, $bond_length_c_up ] );

    # Places third atom on z-axis by the length of the bond in the opposite
    # direction.
    my $bond_length_c_down =
	$ATOMS{$atom_names[0]}{"covalent_radius"}{"length"}[0]
      + $ATOMS{$atom_names[1]}{"covalent_radius"}{"length"}[0];
    push( @atom_coord, [ 0, 0, - $bond_length_c_down ] );

    return \@atom_coord;
}

1;
