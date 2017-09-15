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
    	{ "label_seq_id" => [ $residue_id ] } ),
	[ "label_atom_id" ] );
    my @atom_names = map { $_->[0] } @{ $atom_names };

    # Iterates through every atom of certain residue name and rotates to
    # specified dihedral angle.
    # my %pseudo_rotamers;
    my %generated_rotamers = %{ $atom_site };

    for my $atom_name ( @atom_names ) {
    	# Defines dihedral angles by %ROTATABLE_BONDS description and specified
    	# angles in $angles.
    	my %angles;
    	if( defined rotatable_bonds()->{"$residue_name"}{"$atom_name"} ) {
    	    my $atom_id =
    	    	select_atom_data(
    		filter_atoms( $atom_site,
    		{ "label_atom_id" => [ $atom_name ] } ),
    	    	[ "id" ] )->[0][0];
    	    for my $angle_id
    		( 0..scalar( @{ rotatable_bonds()->{"$residue_name"}
    				                   {"$atom_name"} } ) - 2 ) {
    		    $angles{"chi$angle_id"} =
    			[ $angle_values->{"$residue_id"}{"chi$angle_id"} ];
    	    }

    	    %generated_rotamers =
    		( %generated_rotamers,
    		  %{ generate_pseudo( \%generated_rotamers,
    				      { "id" => [ $atom_id ] },
    				      \%angles ) } );
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
#     $small_angle - angle by which rotation is made.
#     $conformations - possible sidechain movements described by sidechain modeling
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
    my $conformations = $args->{"conformations"};
    my $interactions = $args->{"interactions"};

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
    	# rotatable_bonds.csv parameter file.
    	my @sorted_names =
    	    sort{ scalar( @{ rotatable_bonds->{"$residue_name"}{$a} } )
    	      cmp scalar( @{ rotatable_bonds->{"$residue_name"}{$b} } ) }
    	    keys %{ rotatable_bonds->{"$residue_name"} };

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
		scalar( @{ rotatable_bonds->{"$residue_name"}
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
		print Dumper $angles;
		my %angles =
		    map { ( "chi$_" => [ $angles->[$_], $angles->[$_] ] ) } 0..$angle_count;

		my $pseudo_atom_site =
		    generate_pseudo( $atom_site,
				     { "id" => [ "$current_atom_id" ] },
				     \%angles );
		my $pseudo_atom_id =
		    select_atom_data( $pseudo_atom_site, [ "id" ] )->[0][0];
		print $pseudo_atom_id;
		# LAST:
		# $interactions->( { %{ $atom_site }, %{ $pseudo_atom_site } },
		# 		 { "id" => [ $pseudo_atom_id ] } );
		print Dumper $pseudo_atom_site;
		if( ! exists $pseudo_atom_site->{"$pseudo_atom_id"}{"clashes"} ){
		    push( @allowed_angles, $angles );
		} else {
		    # print Dumper $pseudo_atom_site;
		}
	    }
	    print Dumper \@allowed_angles;
	    die "No possible rotamer solutions were detected."
		if scalar( @allowed_angles ) == 0;
	}
    }

    # return $atom_site;
}

1;
