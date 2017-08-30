package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( generate_pseudo generate_rotamer);

use List::Util qw( max );

use lib qw( ./ );
use PDBxParser qw( filter_atoms select_atom_data );
use Combinatorics qw( permutation );
use LinearAlgebra qw( evaluate_matrix matrix_product );
use LoadParams qw( rotatable_bonds );
use Measure qw( all_dihedral );
use SidechainModels qw( rotation_only );

use Data::Dumper;

# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates pseudo-atoms from side chain models that are written in equation
# form.
# Input:
#     $atom_site - atom site data structure (see PDBxParser).
#     $atom_specifier - hash of hashes for selecting atoms by attributes (see
#     PDBxParser).
#     $angle_values - hash of arrays that describe possible values of dihedral
#     angles.
#     Ex.: { "chi0" => [ 0, pi, 1.5 * $pi, 2 * $pi ],
#            "chi1" => [ 0, 2 * $pi ] }
# Output:
#     $atom_site - atom site data structure with additional pseudo-atoms.
#

sub generate_pseudo
{
    my ( $atom_site, $atom_specifier, $angle_values ) = @_;

    # Determines last id from set of atoms. It will be used for attaching id's
    # to generated pseudo-atoms.
    my $atom_ids = select_atom_data( $atom_site, [ "id" ] );
    my $last_atom_id = max( sort { $a <=> $b } map { $_->[0] } @{ $atom_ids } );

    # Generates model for selected atoms.
    $atom_site = rotation_only( $atom_site );

    my @target_atom_ids =
	map { $_->[0]->[0] }
	select_atom_data(
	filter_atoms( $atom_site, $atom_specifier ),
	[ "id" ] );

    my @angle_values;  # Will be used for temporarily storing sets of angles that
                       # angle permutations will be generated from.
    my @angle_names;
    my $current_resi_id;
    my %current_angles;

    for my $atom_id ( @target_atom_ids ) {
    	@angle_names = sort( { $a cmp $b } keys %{ $angle_values } );
    	undef @angle_values;

    	# Calculates current dihedral angles of rotatable bonds. Will be used
    	# for reseting dihedral angles to 0 degree angle.
	$current_resi_id = $atom_site->{"$atom_id"}{"label_seq_id"};
    	%current_angles =
    	    %{ all_dihedral(
	       filter_atoms( $atom_site,
			     { "label_seq_id" => [ $current_resi_id ] } ) ) };

    	# Iterates through combinations of angles and evaluates conformational
    	# model.
    	for my $angle_name ( @angle_names ) {
    	    push( @angle_values,
    		  [ map
    		    { $_ - $current_angles{"$current_resi_id"}{"$angle_name"} }
    		    @{ $angle_values->{"$angle_name"} } ] );
    	}

    	my %angle_values;
    	my $conf_model = $atom_site->{"$atom_id"}{"conformation"};
    	my $transf_atom_coord;

    	for my $angle_comb (
    	    @{ permutation( scalar( @angle_names ), [], \@angle_values, [] ) } ){
    	    %angle_values =
    		map { $angle_names[$_], $angle_comb->[$_] } 0..$#angle_names;
	    # Converts matrices to GiNaC compatable format and evaluates them.
    	    $transf_atom_coord =
    	    	evaluate_matrix( matrix_product( $conf_model ), \%angle_values );
    	    # Adds generated pseudo-atom to $atom_site.
    	    $last_atom_id++;
    	    %{ $atom_site->{$last_atom_id} } = %{ $atom_site->{$atom_id} };
    	    # Overwrites atom id.
    	    $atom_site->{$last_atom_id}{"id"} = $last_atom_id;
    	    # Overwrites exsisting coordinate values.
    	    $atom_site->{$last_atom_id}{"Cartn_x"} = $transf_atom_coord->[0][0];
    	    $atom_site->{$last_atom_id}{"Cartn_y"} = $transf_atom_coord->[1][0];
    	    $atom_site->{$last_atom_id}{"Cartn_z"} = $transf_atom_coord->[2][0];
    	    # Adds information about used dihedral angles.
    	    $atom_site->{$last_atom_id}{"dihedral_angles"} = \%angle_values;
    	    # Adds additional pseudo-atom flag for future filtering.
    	    $atom_site->{$last_atom_id}{"is_pseudo_atom"} = 1;
    	}
    }

    return $atom_site;
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

    	    # TODO: remove redundant rotation_only function. Should be only
    	    # inside generate_pseudo function. Also, should try to work on
    	    # code readability.
	    # %generated_rotamers =
	    # 	%{ generate_pseudo( rotation_only( \%generated_rotamers ), { "id" => [ $atom_id ] }, \%current_angles ) };
	    print Dumper \%current_angles;
    	    # %generated_rotamers =
    	    # 	%{ generate_pseudo( rotation_only( \%generated_rotamers ), { "id" => [ $atom_id ] }, \%current_angles );
    	}
    }

    # return \%generated_rotamers;
}

1;
