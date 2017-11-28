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
                bond_angle
                bond_length );
use MoleculeProperties qw( %ROTATABLE_BONDS
                           %HYBRIDIZATION );
use PDBxParser qw( create_pdbx_entry
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
    		map { ( $angle_names[$_] => $angle_comb->[$_] ) } 0..$#angle_names;
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

    for my $atom_id ( keys %atom_site ) {
	my $residue_name = $atom_site{$atom_id}{"label_comp_id"};
	my $atom_type = $atom_site{$atom_id}{"type_symbol"};
	my $atom_name = $atom_site{$atom_id}{"label_atom_id"};

	my @connection_ids = @{ $atom_site{"$atom_id"}{"connections"} };
	my @connection_atoms =
	    map { $atom_site{"$_"}{"type_symbol"} } @connection_ids;

	my $hybridization = $HYBRIDIZATION{$residue_name}{$atom_name};
	my $lone_pair_count = $ATOMS{$atom_type}{"lone_pairs"};

	my $hydrogen_count = 4 - scalar( @connection_ids ) - $lone_pair_count;

	# Skips hydrogen addition if there are no hydrogens to add.
	if( $hydrogen_count == 0 ) { next; };

	# Generates name according to the name of heavy atom.
	my @hydrogen_names;
	for my $hydrogen_id ( 1..$hydrogen_count ) {
	    ( my $hydrogen_name = $atom_name ) =~
		s/$atom_type(.?)/H$1$hydrogen_id/;
	    push( @hydrogen_names, $hydrogen_name );
	}

	#            sp3                       sp2               sp
	#
	#            Up(2)                     Up(2)	         Up(2)
	# z          |		               |	         |
	# |_y        C(1) __ Right(3)          C(1)	         C(1)
	# /         / \		              / \	         |
	# x    Left(4) Back(5)	         Left(4) Right(3)        Down(3)
	#
    	# Depending on hybridization and present bond connections, adds missing
    	# hydrogens.
    	# TODO: in the future, should adjust sp3, sp2 angles according to
    	# experimental data, not model.
	# Up atom is a heavy atom that current atom is connected to.
    	if( $hybridization eq "sp3" ) {
    	    if( $hydrogen_count == 3 ) {
    		# Calculates length of bonds.
    		my $right_bond_length =
    		    $ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    		  + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];
    		my $left_bond_length =
    		    $ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    	          + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];
    		my $back_bond_length =
    		    $ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    		  + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

    		# Creates coordinates for hydrogen atoms where the center point
    		# is (0, 0, 0).
    		my @right_atom_coord =
    		    ( [ $right_bond_length * sin( 109.5 * pi() / 180 ) ],
    		      [ 0 ],
    		      [ $right_bond_length * cos( 109.5 * pi() / 180 ) ],
    		      [ 1 ] );
    		my @left_atom_coord =
    		    ( [ $left_bond_length
    		      * cos( 120 * pi() / 180 )
    		      * sin( 109.5 * pi() / 180 ) ],
    		      [ $left_bond_length
    		      * sin( 120 * pi() / 180 )
    		      * sin( 109.5 * pi() / 180 ) ],
    		      [ $left_bond_length
    		      * cos( 109.5 * pi() / 180 ) ],
    		      [ 1 ] );
    		my @back_atom_coord =
    		    ( [ $back_bond_length
    		      * cos( 240 * pi() / 180 )
    		      * sin( 109.5 * pi() / 180 ) ],
    		      [ $back_bond_length
    		      * sin( 240 * pi() / 180 )
    		      * sin( 109.5 * pi() / 180 ) ],
    		      [ $back_bond_length
    		      * cos( 109.5 * pi() / 180 ) ],
    		      [ 1 ] );
    	    } elsif( $hydrogen_count == 2 ) {
		# Decreases bond angle, if lone pairs are present.
		my $right_bond_angle;
		if( $lone_pair_count > 0 ) {
		    $right_bond_angle =
			( 109.5 - $lone_pair_count * 2.5 ) * pi() / 180;
		} else {
		    $right_bond_angle =
			bond_angle(
			    [ [ $atom_site{$connection_ids[0]}{"Cartn_x"},
				$atom_site{$connection_ids[0]}{"Cartn_y"},
				$atom_site{$connection_ids[0]}{"Cartn_z"} ],
			      [ $atom_site{$atom_id}{"Cartn_x"},
				$atom_site{$atom_id}{"Cartn_y"},
				$atom_site{$atom_id}{"Cartn_z"} ],
			      [ $atom_site{$connection_ids[1]}{"Cartn_x"},
				$atom_site{$connection_ids[1]}{"Cartn_y"},
				$atom_site{$connection_ids[1]}{"Cartn_z"} ] ] );
		}
    	    # 	my @left_atom_coord =
    	    # 	    ( [ $left_bond_length
    	    # 	      * cos( 120 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 ) ],
    	    # 	      [ $left_bond_length
    	    # 	      * sin( 120 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 )],
    	    # 	      [ $left_bond_length
    	    # 	      * cos( 109.5 * pi() / 180 )],
    	    # 	      [ 1 ] );
    	    # 	my @back_atom_coord =
    	    # 	    ( [ $back_bond_length
    	    # 	      * cos( 240 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 ) ],
    	    # 	      [ $back_bond_length
    	    # 	      * sin( 240 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 ) ],
    	    # 	      [ $back_bond_length
    	    # 	      * cos( 109.5 * pi() / 180 ) ],
    	    # 	      [ 1 ] );
	    } elsif( $hydrogen_count == 1 ) {
    	    # 	my @back_atom_coord =
    	    # 	    ( [ $back_bond_length
    	    # 	      * cos( 240 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 ) ],
    	    # 	      [ $back_bond_length
    	    # 	      * sin( 240 * pi() / 180 )
    	    # 	      * sin( 109.5 * pi() / 180 ) ],
    	    # 	      [ $back_bond_length
    	    # 	      * cos( 109.5 * pi() / 180 ) ],
    	    # 	      [ 1 ] );
	    }

    	} elsif( $hybridization eq "sp2" ) {
    	    if( $hydrogen_count == 2 ) {
    		# Calculates length of bonds.
    		my $right_bond_length =
    		    $ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    		  + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];
    		my $left_bond_length =
    		    $ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    	          + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

    		# Creates coordinates for hydrogen atoms where the center point
    		# is (0, 0, 0).
    		my @right_atom_coord =
    		    ( [ $right_bond_length * sin( 120 * pi() / 180 ) ],
    		      [ 0 ],
    		      [ $right_bond_length * cos( 120 * pi() / 180 ) ],
    		      [ 1 ] );
    		my @left_atom_coord =
    		    ( [ $left_bond_length * sin( 240 * pi() / 180 ) ],
    		      [ 0 ],
    		      [ $left_bond_length * cos( 240 * pi() / 180 ) ],
    		      [ 1 ] );
    	    }

    	} elsif( $hybridization eq "sp" ) {
    	    # Calculates length of bonds.
    	    my $down_bond_length =
    		$ATOMS{$atom_name}{"covalent_radius"}{"length"}[0]
    	      + $ATOMS{"H"}{"covalent_radius"}{"length"}[0];

    	    # Creates coordinates for hydrogen atoms where the center point
    	    # is (0, 0, 0).
    	    my @down_atom_coord =
    		( [ 0 ],
    		  [ 0 ],
    		  [ - $down_bond_length ],
    		  [ 1 ] );
    	}
    }

    return \%hydrogen_site;
}

1;
