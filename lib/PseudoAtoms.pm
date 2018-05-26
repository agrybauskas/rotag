package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( add_hydrogens
                     generate_library
                     generate_pseudo
                     generate_rotamer
                     library_to_csv );

use List::Util qw( max );
use List::MoreUtils qw( any
                        uniq );
use Math::Trig qw( acos );

use AtomInteractions qw( hard_sphere
                         soft_sphere
                         exponential
                         leonard_jones
                         combined );
use AtomProperties qw( %ATOMS
                       %HYDROGEN_NAMES );
use Combinatorics qw( permutation );
use ConnectAtoms qw( connect_atoms
                     grid_box
                     is_neighbour
                     is_second_neighbour );
use LinearAlgebra qw( find_euler_angles
                      mult_matrix_product
                      pi
                      switch_ref_frame );
use Measure qw( all_dihedral
                dihedral_angle
                bond_angle
                bond_length );
use MoleculeProperties qw( hybridization
                           rotatable_bonds );
use PDBxParser qw( create_pdbx_entry
                   filter
                   to_pdbx );
use Sampling qw( sample_angles );
use SidechainModels qw( rotation_only );

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
    my ( $atom_site, $atom_specifier, $angle_values, $last_atom_id ) = @_;

    $last_atom_id //= max( keys %{ $atom_site } );

    my %atom_site = %{ $atom_site }; # Copy of $atom_site.
    my %pseudo_atom_site;

    my @atom_ids =
	@{ filter( { "atom_site" => \%atom_site,
		     "include" => $atom_specifier,
		     "data" => [ "id" ],
		     "is_list" => 1 } ) };

    for my $atom_id ( @atom_ids ) {
    	# Calculates current dihedral angles of rotatable bonds. Will be used
    	# for reseting dihedral angles to 0 degree angle.
    	my $residue_id = $atom_site{"$atom_id"}{"label_seq_id"};

    	my %angles =
	    %{ all_dihedral(
		   filter( { "atom_site" => \%atom_site,
			     "include" => { "label_seq_id" => [ $residue_id ],
					    "label_alt_id" => [ "." ] } } ) ) };

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

    	my $conformation = $atom_site{"$atom_id"}{"conformation"};
    	for my $angle_comb ( # Abreviation of angle combinations.
    	    @{ permutation( scalar( @angle_names ), [], \@angle_values, [] ) } ){
    	    my %angle_values =
    		map { ( $angle_names[$_] => $angle_comb->[$_] ) }
	        0..$#angle_names;
    	    # Evaluates matrices.
    	    my ( $transf_atom_coord ) =
    	    	@{ mult_matrix_product( $conformation, \%angle_values ) };
    	    # Adds necessary PDBx entries to pseudo atom site.
    	    $last_atom_id++;
    	    create_pdbx_entry(
    	    	{ "atom_site" => \%pseudo_atom_site,
    	    	  "id" => $last_atom_id,
    	    	  "type_symbol" => $atom_site{$atom_id}{"type_symbol"},
    	    	  "label_atom_id" => $atom_site{$atom_id}{"label_atom_id"},
    	    	  "label_alt_id" => "1",
    	    	  "label_comp_id" => $atom_site{$atom_id}{"label_comp_id"},
    	    	  "label_asym_id" => $atom_site{$atom_id}{"label_asym_id"},
    	    	  "label_entity_id" => $atom_site{$atom_id}{"label_entity_id"},
    	    	  "label_seq_id" => $residue_id,
    	    	  "cartn_x" => sprintf( "%.3f", $transf_atom_coord->[0][0] ),
    	    	  "cartn_y" => sprintf( "%.3f", $transf_atom_coord->[1][0] ),
    	    	  "cartn_z" => sprintf( "%.3f", $transf_atom_coord->[2][0] ),
    		  "auth_seq_id" => $atom_site{$atom_id}{"auth_seq_id"} } );
    	    # Adds information about used dihedral angle values and names.
	    $pseudo_atom_site{$last_atom_id}{"dihedral_names"} = \@angle_names;
    	    $pseudo_atom_site{$last_atom_id}{"dihedral_angles"} =
    	    	{ map { ( $_ => $angle_values{$_}
	    		      + $angles{"$residue_id"}{$_} ) }
	    	  @angle_names };
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
    my ( $atom_site, $angle_values, $last_atom_id ) = @_;

    $last_atom_id //= max( keys %{ $atom_site } );

    my %atom_site = %{ $atom_site };
    my %rotamer_atom_site;

    for my $residue_id ( keys %{ $angle_values } ) {
    	my $residue_site =
    	    filter( { "atom_site" => \%atom_site,
    		      "include" => { "label_seq_id" => [ $residue_id ] } } );

    	my $rotatable_bonds = rotatable_bonds( \%atom_site );

    	for my $atom_id ( sort { $a <=> $b } keys %{ $residue_site } ) {
    	    if( ! exists $rotatable_bonds->{$atom_id} ) { next; }

    	    my %angles;
    	    for my $angle_name ( keys %{ $rotatable_bonds->{$atom_id} } ) {
    	    	$angles{$angle_name} =
		    [ $angle_values->{"$residue_id"}{$angle_name} ];
    	    }

    	    %rotamer_atom_site =
    	    	( %rotamer_atom_site,
    	    	  %{ generate_pseudo( { %atom_site, %rotamer_atom_site },
    	    			      { "id" => [ $atom_id ] },
    	    			      \%angles,
    	    			      $last_atom_id ) } );
    	    $last_atom_id++;
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
    my $conf_model = $args->{"conf_model"};
    my $interactions = $args->{"interactions"};
    my $energy_cutoff_atom = $args->{"energy_cutoff_atom"};
    my $energy_cutoff_summed = $args->{"energy_cutoff_summed"};

    $energy_cutoff_summed //= "Inf";

    my %atom_site = %{ $atom_site }; # Copy of $atom_site.

    my %rotamer_library;

    # Generates conformational models before checking for clashes/interactions
    # for given residues.
    if( $conf_model eq "rotation_only" ) {
	rotation_only( filter( { "atom_site" => \%atom_site,
				 "include" =>
				     { "label_seq_id" => $residue_ids } } ) );
    } else {
	die "Conformational model was not defined.";
    }

    # Selection of potential function.
    my $potential_function;
    $potential_function = \&hard_sphere   if $interactions eq "hard_sphere";
    $potential_function = \&soft_sphere   if $interactions eq "soft_sphere";
    $potential_function = \&exponential   if $interactions eq "exponential";
    $potential_function = \&leonard_jones if $interactions eq "leonard_jones";
    $potential_function = \&combined      if $interactions eq "combined";

    # Creates the grid box that has edge length of sum of all bonds of the
    # longest side-chain branch in arginine. Length: 3 * (C-C) + (C-N) + 2
    # * (C=N) + (N-H).
    # TODO: should consider using shorter distances, because bonds have limits
    # on maximum bending and having shorter edge length reduces calculation time.
    my ( $grid_box, undef ) =
	grid_box( \%atom_site,
		  7 * $ATOMS{"C"}{"covalent_radius"}{"length"}->[0]
		+ 2 * $ATOMS{"N"}{"covalent_radius"}{"length"}->[0]
	        + 3 * $ATOMS{"C"}{"covalent_radius"}{"length"}->[1]
		+ 3 * $ATOMS{"N"}{"covalent_radius"}{"length"}->[1]
		+     $ATOMS{"H"}{"covalent_radius"}{"length"}->[0] );

    # Finds where CA of target residues are.
    my %target_cell_idxs;
    for my $cell_idx ( keys %{ $grid_box } ) {
	for my $atom_id ( @{ $grid_box->{"$cell_idx"} } ) {
	    my $residue_id = $atom_site->{$atom_id}{"label_seq_id"};
	    my $atom_name = $atom_site->{$atom_id}{"label_atom_id"};
	    if( $atom_name eq 'CA'
	     && any { $residue_id eq $_ } @{ $residue_ids } ) {
		push( @{ $target_cell_idxs{$cell_idx} }, $residue_id );
	    }
	}
    }

    for my $cell_idxs ( sort { $a cmp $b } keys %target_cell_idxs ) {
    	my @cell_idxs = split( ",", $cell_idxs );
    	my @neighbour_atom_ids; # The array will contain all atoms of the
                                # neighbouring 26 cells.

    	# $i represents x, $j - y, $k - z coordinates.
    	for my $i ( ( $cell_idxs[0] - 1..$cell_idxs[0] + 1 ) ) {
    	for my $j ( ( $cell_idxs[1] - 1..$cell_idxs[1] + 1 ) ) {
    	for my $k ( ( $cell_idxs[2] - 1..$cell_idxs[2] + 1 ) ) {
    	if( exists $grid_box->{"$i,$j,$k"} ) {
    	    push( @neighbour_atom_ids, @{ $grid_box->{"$i,$j,$k"} } ); } } } }

    	for my $residue_id ( @{ $target_cell_idxs{$cell_idxs} } ) {
    	    my $residue_site =
    		filter( { "atom_site" => \%atom_site,
    			  "include" => { "label_seq_id" => [ $residue_id ] } } );

    	    my $rotatable_bonds = rotatable_bonds( $residue_site );
    	    if( ! %{ $rotatable_bonds } ) { next; }

    	    # Because the change of side-chain position might impact the
    	    # surrounding, iteraction site consists of only main chain atoms.
	    # TODO: decide, if interactions should be controlled by command line.
    	    my %interaction_site =
    		%{ filter( { "atom_site" => \%atom_site,
    			     "include" =>
    			         { "id" => \@neighbour_atom_ids,
    				   "label_atom_id" => [ "N", "CA", "C", "O",
    							"OXT", "CB", "H", "H2",
    							"HA2", "HA3", "HB1",
    							"HB2", "HB3", "HXT"
    				                      ] } } ) };

    	    # Goes through each atom in side chain and calculates interaction
    	    # potential with surrounding atoms. CA and CB are non-movable atoms
    	    # so, they are marked as starting atoms.
    	    my $ca_atom_id =
    	    	filter( { "atom_site" => $residue_site,
    	    		  "include" => { "label_atom_id" => [ "CA" ] },
    	    		  "data" => [ "id" ] } )->[0][0];
    	    my $cb_atom_id =
    	    	filter( { "atom_site" => $residue_site,
    	    		  "include" => { "label_atom_id" => [ "CB" ] },
    	    		  "data" => [ "id" ] } )->[0][0];
    	    my @visited_atom_ids = ( $ca_atom_id, $cb_atom_id );
    	    my @next_atom_ids =
    	    	grep { $_ ne $ca_atom_id }
    	    	@{ $residue_site->{$cb_atom_id}{"connections"} };

    	    my @sampled_angles =
    	    	map { [ $_ ] }
    	        @{ sample_angles( [ [ 0, 2 * pi() ] ], $small_angle ) };
    	    my @allowed_angles = @sampled_angles;

            # @zero_energies is a helper variable for permutation() in order to
            # mimick the permutated angles and match their values correctly.
    	    my @zero_energies = map { [ 0 ] } @sampled_angles;
            my @allowed_energies = @zero_energies;

    	    while( scalar( @next_atom_ids ) != 0 ) {
    	    	my @neighbour_atom_ids;
    	    	for my $atom_id ( @next_atom_ids ) {
    	    	    # Adds more angle combinations if there are more than one
    	    	    # rotatable bonds.
    	    	    if( scalar( @{ $allowed_angles[0] } )
    	    	      < scalar( keys %{ $rotatable_bonds->{$atom_id} } ) ) {
    	    	    	@allowed_angles =
    	    	    	    @{ permutation( 2, [], [ \@allowed_angles,
    	    	    				     \@sampled_angles ], [] ) };
    	    	    	@allowed_energies =
    	    	    	    @{ permutation( 2, [], [ \@allowed_energies,
    	    	    				     \@zero_energies ], [] ) };
    	    		# Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
    	    		@allowed_angles =
    	    		    map { [ @{ $_->[0] }, @{ $_->[1] } ] }
    	    		    @allowed_angles;
    	    		@allowed_energies =
    	    		    map { [ $_->[0][0] ] }
    	    		    @allowed_energies;
    	    	    }

    	    	    # Marks visited atoms.
    	    	    push( @visited_atom_ids, $atom_id );

    	    	    # Marks neighbouring atoms.
    	    	    push( @neighbour_atom_ids,
    	    		  @{ $atom_site{$atom_id}{"connections"} } );

    	    	    # Starts calculating potential energy.
    	    	    my @next_allowed_angles;
    	    	    my @next_allowed_energies;
                    for( my $i = 0; $i <= $#allowed_angles; $i++ ) {
                        my $angles = $allowed_angles[$i];
                        my $energies = $allowed_energies[$i]->[0];
    	    		my %angles =
    	    		    map { ( "chi$_" => [ $angles->[$_] ] ) }
    	    		    0..$#{ $angles };
    	    		my $pseudo_atom_site =
    	    		    generate_pseudo( \%atom_site,
    	    				     { "id" => [ "$atom_id" ] },
    	    				     \%angles );
    	    		my $pseudo_atom_id = ( keys %{ $pseudo_atom_site } )[0];
    	    		my $pseudo_origin_id =
    	    		    $pseudo_atom_site->{$pseudo_atom_id}
                                               {"origin_atom_id"};

			# TODO: decide, if choosing that each interaction should
			# not reach certain boundaries will be suitable for the
			# model. Otherwise, there might be interactions that
			# outweight the clashes and depending on the order of
			# calculations might produce different results. Or maybe,
			# the sum of potential should be calculated.
			my $potential_energy;
                        my $potential_sum;
    	    		foreach my $interaction_id ( keys %interaction_site ) {
    	    		    if( ( ! is_neighbour( \%atom_site,
    	    					  $pseudo_origin_id,
    	    					  $interaction_id ) )
    	    		     && ( ! is_second_neighbour( \%atom_site,
    	    		     				 $pseudo_origin_id,
    	    		     				 $interaction_id ) )
    	    			) {
    	    			$potential_energy =
    	    			    $potential_function->(
    	    				$pseudo_atom_site->{$pseudo_atom_id},
    	    				$atom_site{$interaction_id} );
                                $potential_sum += $potential_energy;
    	    			last if $potential_energy > $energy_cutoff_atom;
    	    		    }
    	    		}

    	    		# Writes allowed angles to @next_allowed_angles that will
    	    		# be passed to more global @allowed_angles. Checks the
			# last calculated potential. If potential was greater
			# than the cutoff, then calculation was halted, but the
			# value remained.
    	    		if( $potential_energy <= $energy_cutoff_atom ) {
    	    		    push( @next_allowed_angles, $angles );
    	    		    push( @next_allowed_energies,
                                  [ $energies + $potential_sum ] );
    	    		}
    	    	    }

    	    	    if( scalar( @allowed_angles ) > 0 ) {
    	    		@allowed_angles = @next_allowed_angles;
    	    		@allowed_energies = @next_allowed_energies;
    	    	    } else {
    	    		die "No possible rotamer solutions were detected.";
    	    	    }
    	    	}

    	    	# Determines next atoms that should be visited.
    	    	@next_atom_ids = (); # Resets value for the new ones to be
              	                     # appended.
    	    	for my $neighbour_atom_id ( uniq @neighbour_atom_ids ) {
    	    	    if( ( ! grep { $neighbour_atom_id eq $_ }
    	    		    @visited_atom_ids ) ) {
    	    		push( @next_atom_ids, $neighbour_atom_id );
    	    	    }
    	    	}
    	    }

    	    # TODO: remember to add check on inter-atom clashing inside
    	    # side-chain itself.
            for( my $i = 0; $i <= $#allowed_angles; $i++ ) {
    	    	my $angles = $allowed_angles[$i];
                my $energies = $allowed_energies[$i]->[0];
                if( $energies <= $energy_cutoff_summed ) {
                    push( @{ $rotamer_library{"$residue_id"} },
                          { "angles" => { map { ( "chi$_" => $angles->[$_] ) }
                                              ( 0..$#{ $angles } ) },
                            "potential" => $interactions,
                            "potential_energy_value" => $energies } );
                }
    	    }
    	}
    }

    return \%rotamer_library;
}

sub add_hydrogens
{
    my ( $atom_site ) = @_;

    my %atom_site = %{ $atom_site };

    connect_atoms( \%atom_site );
    hybridization( \%atom_site );

    my %hydrogen_site;
    my $last_atom_id = max( keys %{ $atom_site } );

    for my $atom_id ( sort { $a <=> $b } keys %atom_site ) {
    	my $atom_type = $atom_site{$atom_id}{"type_symbol"};
    	my $atom_name = $atom_site{$atom_id}{"label_atom_id"};

	my $residue_name = $atom_site{$atom_id}{"label_comp_id"};
	my $residue_id = $atom_site{$atom_id}{"label_seq_id"};
	my $residue_site =
	    filter( { "atom_site" => $atom_site,
		      "include" => { "label_seq_id" => [ $residue_id ] } } );
	my %residue_coord =
	    %{ filter( { "atom_site" => $residue_site,
			 "data" => [ "Cartn_x", "Cartn_y", "Cartn_z" ],
			 "data_with_id" => 1 } ) };

    	my $hydrogen_names = $HYDROGEN_NAMES{$residue_name}{$atom_name};

    	if( ! $hydrogen_names ) { next; }; # Exits early if there should be no
	                                   # hydrogens connected to the atom.

    	my @connection_ids = @{ $atom_site{"$atom_id"}{"connections"} };
    	my @connection_names =
    	    map { $atom_site{"$_"}{"label_atom_id"} } @connection_ids;

    	my $hybridization = $residue_site->{$atom_id}{"hybridization"};
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
			$residue_site->{$second_neighbour}{"hybridization"};
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

# ------------------------------- Library to STDOUT --------------------------- #

sub library_to_csv
{
    my ( $rotamer_library, $angle_units ) = @_;

    $angle_units //= "radians";

    print( "\"rotamer_id\",\"residue_id\",\"angle_name\",\"angle_value\"\n" );

    my $rotamer_id = 1;
    for my $residue_id ( keys %{ $rotamer_library } ) {
	for my $rotamer ( @{ $rotamer_library->{$residue_id} } ) {
	    for my $angle_name ( sort { $a cmp $b } keys %{ $rotamer } ) {
		if( $angle_units eq "radians" ) {
		    printf( "%d,%d,\"%s\",%.3f\n",
			    $rotamer_id,
			    $residue_id,
			    $angle_name,
			    $rotamer->{$angle_name} );
		} elsif( $angle_units eq "degrees" ) {
		    printf( "%d,%d,\"%s\",%.2f\n",
			    $rotamer_id,
			    $residue_id,
			    $angle_name,
			    $rotamer->{$angle_name} * 180 / pi() );
		}
	    }
	    $rotamer_id++;
	}
    }

    return;
}

1;
