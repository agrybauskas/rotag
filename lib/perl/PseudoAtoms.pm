package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( generate_library generate_pseudo generate_rotamer );

use List::Util qw( max );

use lib qw( ./ );
use PDBxParser qw( filter_atoms select_atom_data );
use Combinatorics qw( permutation );
use LinearAlgebra qw( evaluate_matrix matrix_product pi );
use LoadParams qw( rotatable_bonds );
use Measure qw( all_dihedral );
use Sampling qw( sample_angles );
use SidechainModels qw( rotation_only );
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
			     { "label_seq_id" => [ $residue_id ] } ) ) };

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
    	    	evaluate_matrix( matrix_product( $conformation ),
				 \%angle_values );
    	    # TODO: decide what data should be copied and what data should be
    	    # with ? or . symbols in $atom_site for pseudo-atoms.
    	    # Adds generated pseudo-atom to $atom_site.
    	    $last_atom_id++;
    	    $pseudo_atom_site{$last_atom_id}{"id"} = $last_atom_id;
    	    # Adds atom type.
    	    $pseudo_atom_site{$last_atom_id}{"type_symbol"} =
		$atom_site->{$atom_id}{"type_symbol"};
    	    # Adds coordinate values.
    	    $pseudo_atom_site{$last_atom_id}{"Cartn_x"} =
		$transf_atom_coord->[0][0];
    	    $pseudo_atom_site{$last_atom_id}{"Cartn_y"} =
		$transf_atom_coord->[1][0];
    	    $pseudo_atom_site{$last_atom_id}{"Cartn_z"} =
		$transf_atom_coord->[2][0];
    	    # Adds information about used dihedral angles.
    	    $pseudo_atom_site{$last_atom_id}{"dihedral_angles"} =
		\%angle_values;
    	    # Adds additional pseudo-atom flag for future filtering.
    	    $pseudo_atom_site{$last_atom_id}{"is_pseudo_atom"} = 1;
    	}
    }

    return \%pseudo_atom_site;
}

#
# Generates rotamers according to given angle values.
# Input:
#     $atom_site - atom site data structure (see PDBxParser).
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
    my @resi_ids = keys %{ $angle_values };
    my $resi_id = $resi_ids[0];
    my $resi_name =
    	select_atom_data(
	    filter_atoms( $atom_site,
	    { "label_seq_id" => [ $resi_id ] } ), [ "label_comp_id" ] )->[0][0];

    my $atom_labels =
    	select_atom_data(
	    filter_atoms( $atom_site,
	    { "label_seq_id" => [ $resi_id ] } ), [ "label_atom_id" ] );
    my @atom_labels = map { $_->[0] } @{ $atom_labels };

    # Iterates through every atom of certain residue name and rotates to
    # specified dihedral angle.
    my %generated_rotamers = %{ $atom_site };
    my $atom_id;
    my %current_angles;

    for my $atom_label ( @atom_labels ) {
    	# Defines dihedral angles by %ROTATABLE_BONDS description and specified
    	# angles in $angles.
    	undef %current_angles;
    	if( defined rotatable_bonds()->{"$resi_name"}{"$atom_label"} ) {
    	    $atom_id =
    	    	select_atom_data(
		filter_atoms( $atom_site,
		{ "label_atom_id" => [ $atom_label ] } ),
    	    	[ "id" ] )->[0][0];
    	    for my $angle_id
    		( 0..scalar( @{ rotatable_bonds()->{"$resi_name"}
    				                   {"$atom_label"} } ) - 2 ) {
    		    $current_angles{"chi$angle_id"} =
    			[ $angle_values->{"$resi_id"}{"chi$angle_id"} ];
	    }

	    %generated_rotamers =
		( %generated_rotamers,
		  %{ generate_pseudo( \%generated_rotamers,
				      { "id" => [ $atom_id ] },
				      \%current_angles ) } );
    	}
    }

    return \%generated_rotamers;
}

#
# Generates rotamer libraries by specified arguments that include atom
# movements and interactions between atoms.
# Input:
#     $atom_site - atom site data structure (see PDBxParser).
#     $residue_ids - array of residue ids.
#     $movements - possible sidechain movements described by sidechain modeling
#     functions in SidechainModels.pm.
#     $interactions - interaction models described by functions in
#     AtomInteractions.pm
# Output:
#     %generated_library - atom site data structure with additional data
#

sub generate_library
{
    my ( $args ) = @_;
    my $atom_site = $args->{"atom_site"};
    my $residue_ids = $args->{"residue_ids"};
    my $small_angle = $args->{"small_angle"};
    my $movements = $args->{"movements"};
    my $interactions = $args->{"interactions"};

    my %generated_library = %{ $atom_site };

    my $resi_site; # Atom site data structure for residue.
    my $resi_name; # Residue name, such as, SER, GLU and etc.
    my $current_atom; # Atom site data structure, but for one atom.
    my $current_atom_id;
    my @sorted_names; # Sorted atom names according to the quantity of rotatable
                      # bonds.

    my $pseudo_id;
    my $pseudo_site;

    my $sampled_angles = sample_angles( [ [ 0, 2 * pi() ] ], $small_angle );
    my @angle_names;
    my %all_current_angles; # TODO: change variable name.
    my @current_angles; # TODO: change variable name.
    my %current_angles; # TODO: change variable name.

    for my $residue_id ( @{ $residue_ids } ) {
	$resi_site =
	    filter_atoms( $atom_site, { "label_seq_id" => [ $residue_id ] } );
	$resi_name =
	    select_atom_data( $resi_site, [ "label_comp_id" ] )->[0][0];

	# Sorts atom names by the quantity of rotatable bonds described in
	# rotatable_bonds.csv parameter file.
	@sorted_names =
	    sort{ scalar( @{ rotatable_bonds->{"$resi_name"}{$a} } )
	      cmp scalar( @{ rotatable_bonds->{"$resi_name"}{$b} } ) }
	    keys %{ rotatable_bonds->{"$resi_name"} };

	# Iterates through sorted atoms and tries to detect interactions.
	undef %current_angles; # Resets angles for each new residue.

	for my $atom_name ( @sorted_names ) {
	    $current_atom =
		filter_atoms( $resi_site,
			      { "label_atom_id" => [ "$atom_name" ] } );
	    $current_atom_id =
	    	select_atom_data( $current_atom, [ "id" ] )->[0][0];
	    @angle_names =
		map { "chi$_" }
	        ( 0..scalar( @{ rotatable_bonds()->{"$resi_name"}
		                    		   {"$atom_name"} } ) - 2 );

	    # Checks for previously checked angles. If not asigned, adds all
	    # possible angles.
	    for my $angle_name ( @angle_names ) {
	    	if( not defined $all_current_angles{"$angle_name"} ) {
	    	    $all_current_angles{"$angle_name"} = $sampled_angles;
	    	}
	    }

	    # Generates combinations of available angles.
	    @current_angles =
	    	map { $all_current_angles{$_} } keys %all_current_angles;

	    for my $angle_comb ( @{ permutation( scalar( @angle_names ), [],
	    					 \@current_angles, [] ) } ) {
	    	%current_angles =
	    	    map { $angle_names[$_], [ $angle_comb->[$_] ] }
	    	    0..$#angle_names;

		# Checks for clashes/interactions. If any detected, removes angle
		# from %all_current_angles.
		# TODO: should be applicable not only to rotation_only model. See
		# generate_pseudo function.
		$pseudo_site =
		    generate_pseudo( \%generated_library,
				     { "id" => [ "$current_atom_id" ] },
				     \%current_angles );
		$pseudo_id = select_atom_data( $pseudo_site, [ "id" ] )->[0][0];

	    }
	}
    }
}

1;
