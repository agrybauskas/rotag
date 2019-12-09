package Measure;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
    our @EXPORT_OK = qw( all_dihedral
                         around_distance
                         bond_angle
                         bond_length
                         dihedral_angle
                         distance
                         distance_squared
                         rmsd
                         rmsd_sidechains
                         energy );
}

use Carp;
use Clone qw( clone );
use Math::Trig;
use List::MoreUtils qw( any
                        uniq );
use List::Util qw( sum );

use AtomProperties qw( sort_atom_names );
use ConnectAtoms qw( is_neighbour
                     is_second_neighbour );
use Grid qw( grid_box
             identify_neighbour_cells );
use Energy;
use ForceField::Bonded;
use ForceField::NonBonded;
use PDBxParser qw( filter
                   filter_by_unique_residue_key
                   split_by );
use LinearAlgebra qw( matrix_sub
                      vector_cross );
use BondProperties qw( rotatable_bonds );
use Version qw( $VERSION );

our $VERSION = $VERSION;

my %potentials = (
    'composite' => {
        'bonded' => {
            'torsion' => \&ForceField::Bonded::torsion_components,
        },
        'non_bonded' => {
            'lennard_jones' => \&ForceField::NonBonded::lennard_jones,
            'coulomb' => \&ForceField::NonBonded::coulomb,
            'h_bond' => \&ForceField::NonBonded::h_bond,
        }
    },
    'torsion' => {
        'bonded' => {
            'torsion' => \&ForceField::Bonded::torsion_components,
        }
    },
    'hard_sphere' => {
        'non_bonded' => {
            'hard_sphere' => \&ForceField::NonBonded::hard_sphere,
        }
    },
    'soft_sphere' => {
        'non_bonded' => {
            'soft_sphere' => \&ForceField::NonBonded::soft_sphere
        }
    },
    'lennard_jones' => {
        'non_bonded' => {
            'lennard_jones' => \&ForceField::NonBonded::lennard_jones,
        }
    },
    'coulomb' => {
        'non_bonded' => {
            'coulomb' => \&ForceField::NonBonded::coulomb,
        }
    },
    'h_bond' => {
        'non_bonded' => {
            'h_bond' => \&ForceField::NonBonded::h_bond,
        }
    }
);

# --------------------------- Molecule parameters ----------------------------- #

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
        sqrt( ( $atom_coord->[1][0] - $atom_coord->[0][0] )**2 +
              ( $atom_coord->[1][1] - $atom_coord->[0][1] )**2 +
              ( $atom_coord->[1][2] - $atom_coord->[0][2] )**2 );

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
                      $atom_coord->[0][2] - $atom_coord->[1][2], );
    my @vector_bc = ( $atom_coord->[2][0] - $atom_coord->[1][0],
                      $atom_coord->[2][1] - $atom_coord->[1][1],
                      $atom_coord->[2][2] - $atom_coord->[1][2], );

    my $vector_product = $vector_ab[0] * $vector_bc[0] +
                         $vector_ab[1] * $vector_bc[1] +
                         $vector_ab[2] * $vector_bc[2];

    my $length_ab = sqrt( $vector_ab[0]**2 +
                          $vector_ab[1]**2 +
                          $vector_ab[2]**2 );
    my $length_bc = sqrt( $vector_bc[0]**2 +
                          $vector_bc[1]**2 +
                          $vector_bc[2]**2 );

    my $bond_angle = acos $vector_product / ( $length_ab * $length_bc );

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
    my $vector_cross_ab = vector_cross( @{ $vector_a }, @{ $vector_b } );
    my $vector_cross_bc = vector_cross( @{ $vector_b }, @{ $vector_c } );

    # Calculates length for each cross product.
    my $vector_length_ab = sqrt $vector_cross_ab->[0]**2 +
                                $vector_cross_ab->[1]**2 +
                                $vector_cross_ab->[2]**2;
    my $vector_length_bc = sqrt $vector_cross_bc->[0]**2 +
                                $vector_cross_bc->[1]**2 +
                                $vector_cross_bc->[2]**2;

    # Calculates normal vectors for each cross product of two vectors.
    my @normal_vector_ab = map { $_ / $vector_length_ab } @{ $vector_cross_ab };
    my @normal_vector_bc = map { $_ / $vector_length_bc } @{ $vector_cross_bc };

    # Finishes orthonormal frame from normal vector ab, vector b and its cross
    # product.
    my $vector_length_b = sqrt $vector_b->[0][0]**2 +
                               $vector_b->[0][1]**2 +
                               $vector_b->[0][2]**2;
    my @normal_vector_b = map { $_ / $vector_length_b } @{ $vector_b->[0] };

    my @orthonormal_cross =
        vector_cross( \@normal_vector_ab, \@normal_vector_b );

    # Using orthonormal frame, projections from vector a and c are
    # generated and angle calculated.
    my $dihedral_angle =
        - atan2 $orthonormal_cross[0][0] * $normal_vector_bc[0] +
                $orthonormal_cross[0][1] * $normal_vector_bc[1] +
                $orthonormal_cross[0][2] * $normal_vector_bc[2],
                $normal_vector_ab[0] * $normal_vector_bc[0] +
                $normal_vector_ab[1] * $normal_vector_bc[1] +
                $normal_vector_ab[2] * $normal_vector_bc[2];

    return $dihedral_angle;
}

#
# Calculates dihedral angles for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it).
# Usage of connect_atoms and hybridization functions are necessary for correct
# calculations.
# Input:
#     $atom_site - atom data structure.
#     $options->{'calc_mainchain'} - additionally calculates phi and psi
#     mainchain dihedral angles.
# Output:
#     $residue_angles - data structure that relates residue id and angle values.
#     Ex.:
#       resi_id, angle_name, angle value
#     { 18 => { 'chi1' => '-3.14' } }
#

sub all_dihedral
{
    my ( $atom_site, $options ) = @_;
    my ( $calc_mainchain, $reference_atom_site ) = (
        $options->{'calc_mainchain'},
        $options->{'reference_atom_site'},
    );

    $calc_mainchain //= 0;
    $reference_atom_site //= $atom_site;

    my %atom_site = %{ $atom_site }; # Copy of $atom_site.

    my $residue_groups =
        split_by( { 'atom_site' => \%atom_site, 'append_dot' => 1 } );

    # Iterates through residue ids and, according to the parameter file,
    # calculates dihedral angles of each rotatable bond.
    my %residue_angles;

    for my $residue_unique_key ( keys %{ $residue_groups } ) {
        my $residue_site =
            filter( { 'atom_site' => \%atom_site,
                      'include' =>
                          { 'id' => $residue_groups->{$residue_unique_key} } } );

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

        my %angle_values;

        # Calculates main-chain phi, psi angles.
        if( $calc_mainchain ) {
            my $n_atom_id =
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'label_atom_id' => [ 'N' ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];
            my $ca_atom_id =
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'label_atom_id' => [ 'CA' ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];
            my $c_atom_id =
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'label_atom_id' => [ 'C' ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];

            # TODO: look if these filter slow down calculations drastically.
            my $prev_c_atom_id;
            if( defined $n_atom_id &&
                defined $residue_site->{$n_atom_id}{'connections'} ) {
                $prev_c_atom_id = filter(
                    { 'atom_site' => $reference_atom_site,
                      'include' =>
                          { 'id' => $residue_site->{$n_atom_id}{'connections'},
                            'label_atom_id' => [ 'C' ] },
                            'data' => [ 'id' ],
                            'is_list' => 1 }
                )->[0];
            }
            my $next_n_atom_id;
            if( defined $c_atom_id &&
                defined $residue_site->{$c_atom_id}{'connections'} ) {
                $next_n_atom_id = filter(
                    { 'atom_site' => $reference_atom_site,
                      'include' =>
                          { 'id' => $residue_site->{$c_atom_id}{'connections'},
                            'label_atom_id' => [ 'N' ] },
                            'data' => [ 'id' ],
                            'is_list' => 1 }
                )->[0];
            }

            # Calculates phi angle if 'C' atom of previous residue is present.
            if( defined $prev_c_atom_id && defined $n_atom_id &&
                defined $ca_atom_id && defined $ca_atom_id ) {
                $angle_values{'phi'}{'atom_ids'} =
                    [ $prev_c_atom_id, $n_atom_id, $ca_atom_id, $c_atom_id ];
                $angle_values{'phi'}{'value'} =
                    dihedral_angle(
                        [ [ $reference_atom_site->{$prev_c_atom_id}{'Cartn_x'},
                            $reference_atom_site->{$prev_c_atom_id}{'Cartn_y'},
                            $reference_atom_site->{$prev_c_atom_id}{'Cartn_z'}, ],
                          [ $atom_site->{$n_atom_id}{'Cartn_x'},
                            $atom_site->{$n_atom_id}{'Cartn_y'},
                            $atom_site->{$n_atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$ca_atom_id}{'Cartn_x'},
                            $atom_site->{$ca_atom_id}{'Cartn_y'},
                            $atom_site->{$ca_atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$c_atom_id}{'Cartn_x'},
                            $atom_site->{$c_atom_id}{'Cartn_y'},
                            $atom_site->{$c_atom_id}{'Cartn_z'} ], ] );
            }

            # Calculates psi angle.
            if( defined $next_n_atom_id && defined $n_atom_id &&
                defined $ca_atom_id && defined $c_atom_id ) {
                $angle_values{'psi'}{'atom_ids'} =
                    [ $n_atom_id, $ca_atom_id, $c_atom_id, $next_n_atom_id ];
                $angle_values{'psi'}{'value'} =
                    dihedral_angle(
                        [ [ $atom_site->{$n_atom_id}{'Cartn_x'},
                            $atom_site->{$n_atom_id}{'Cartn_y'},
                            $atom_site->{$n_atom_id}{'Cartn_z'}, ],
                          [ $atom_site->{$ca_atom_id}{'Cartn_x'},
                            $atom_site->{$ca_atom_id}{'Cartn_y'},
                            $atom_site->{$ca_atom_id}{'Cartn_z'} ],
                          [ $atom_site->{$c_atom_id}{'Cartn_x'},
                            $atom_site->{$c_atom_id}{'Cartn_y'},
                            $atom_site->{$c_atom_id}{'Cartn_z'} ],
                          [ $reference_atom_site->{$next_n_atom_id}{'Cartn_x'},
                            $reference_atom_site->{$next_n_atom_id}{'Cartn_y'},
                            $reference_atom_site->{$next_n_atom_id}{'Cartn_z'} ], ] );
            }
        }

        # Calculates every side-chain dihedral angle.
        for my $angle_name ( keys %uniq_rotatable_bonds ) {
            # First, checks if rotatable bond has fourth atom produce dihedral
            # angle. It is done by looking at atom connections - if rotatable
            # bond ends with terminal atom, then this bond is excluded.
            if( scalar( @{ $residue_site->
                               {$uniq_rotatable_bonds{$angle_name}->[1]}
                               {'connections'} } ) < 2 ){ next; }

            # Chooses proper atom ids for calculating dihedral angles.
            my $second_atom_id = $uniq_rotatable_bonds{$angle_name}->[0];
            my $third_atom_id = $uniq_rotatable_bonds{$angle_name}->[1];
            my @second_connections = # Second atom connections, except third.
                grep { $_ ne $third_atom_id }
                @{ $residue_site->{$second_atom_id}{'connections'} };
            my $first_atom_name =
                sort_atom_names(
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'id' => \@second_connections },
                          'data' => [ 'label_atom_id' ],
                          'is_list' => 1 } ) )->[0];
            my $first_atom_id =
                filter( { 'atom_site' => $residue_site,
                          'include' =>
                              { 'label_atom_id' => [ $first_atom_name ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];
            my @third_connections = # Third atom connections, except second.
                grep { $_ ne $second_atom_id }
                @{ $residue_site->{$third_atom_id}{'connections'} };
            my $fourth_atom_name =
                sort_atom_names(
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'id' => \@third_connections },
                          'data' => [ 'label_atom_id' ],
                          'is_list' => 1 } ) )->[0];
            my $fourth_atom_id =
                filter( { 'atom_site' => $residue_site,
                          'include' =>
                              { 'label_atom_id' => [ $fourth_atom_name ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];

            # Extracts coordinates for dihedral angle calculations.
            my $first_atom_coord =
                [ $residue_site->{$first_atom_id}{'Cartn_x'},
                  $residue_site->{$first_atom_id}{'Cartn_y'},
                  $residue_site->{$first_atom_id}{'Cartn_z'} ];
            my $second_atom_coord =
                [ $residue_site->{$second_atom_id}{'Cartn_x'},
                  $residue_site->{$second_atom_id}{'Cartn_y'},
                  $residue_site->{$second_atom_id}{'Cartn_z'} ];
            my $third_atom_coord =
                [ $residue_site->{$third_atom_id}{'Cartn_x'},
                  $residue_site->{$third_atom_id}{'Cartn_y'},
                  $residue_site->{$third_atom_id}{'Cartn_z'} ];
            my $fourth_atom_coord =
                [ $residue_site->{$fourth_atom_id}{'Cartn_x'},
                  $residue_site->{$fourth_atom_id}{'Cartn_y'},
                  $residue_site->{$fourth_atom_id}{'Cartn_z'} ];

            $angle_values{$angle_name}{'atom_ids'} =
                [ $first_atom_id,
                  $second_atom_id,
                  $third_atom_id,
                  $fourth_atom_id ];
            $angle_values{$angle_name}{'value'} =
                dihedral_angle( [ $first_atom_coord,
                                  $second_atom_coord,
                                  $third_atom_coord,
                                  $fourth_atom_coord ] );
        }

        if( %angle_values ) {
            %{ $residue_angles{$residue_unique_key} } = %angle_values;
        }
    }

    return \%residue_angles;
}

#
# Calculates the squared distance between two atoms.
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm).
# Output:
#     $distance_squared - value of calculated squared distance.
#

sub distance_squared
{
    my ( $atom_i, $atom_j ) = @_;

    my $distance_squared =
        ( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2 +
        ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2 +
        ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    return $distance_squared;
}

#
# Calculates the distance between two atoms.
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm).
# Output:
#     $distance - value of the calculated distance.
#

sub distance
{
    my ( $atom_i, $atom_j ) = @_;

    return sqrt( distance_squared( $atom_i, $atom_j ) );
}

#
# Selects the atoms that are at specified distance from selected atoms.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.obtain_atom_site);
#     $atom_specifier - atom selector data structure (see PDBxParser::filter);
#     $distance - max distance from which atoms should be included.
# Output:
#     %around_atom_site - atom site data structure of selected atoms.
#

sub around_distance
{
    my ( $parameters, $atom_site, $atom_specifier, $distance ) = @_;

    my @atom_ids = @{ filter( { 'atom_site' => $atom_site,
                                'include' => $atom_specifier,
                                'data' => [ 'id' ],
                                'is_list' => 1 } ) };

    # For each cell, checks neighbouring cells. Creates box around atoms, makes
    # grid with edge length of max covalent radii of the parameter file.
    my ( $grid_box, $atom_cell_pos ) =
        grid_box( $parameters, $atom_site, $distance * 2, \@atom_ids );
    my $neighbour_cells = identify_neighbour_cells( $grid_box, $atom_cell_pos );

    # Checks for neighbouring cells for each cell.
    my %around_atom_site;
    foreach my $cell ( keys %{ $atom_cell_pos } ) {
        foreach my $atom_id ( @{ $atom_cell_pos->{$cell} } ) {
            foreach my $neighbour_id ( @{ $neighbour_cells->{$cell} } ) {
        	if( ( ! any { $neighbour_id eq $_ } @atom_ids ) &&
                    ( distance_squared(
                          $atom_site->{$atom_id},
                          $atom_site->{$neighbour_id} ) <= $distance ** 2 ) ) {
        	    $around_atom_site{$neighbour_id} = $atom_site->{$neighbour_id};
        	} }
        }
    }

    return \%around_atom_site;
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

    # Gives the error if the size of the sets are different.
    if( scalar @{ $first_set } != scalar @{ $second_set } ) {
        confess 'comparing different sizes of sets of the atoms is not allowed.';
    }

    # Sums up sqaured differences of coordinates.
    for( my $i = 0; $i <= $#{ $first_set }; $i++ ) {
        $rmsd += ( $first_set->[$i][0] - $second_set->[$i][0] )**2 +
                 ( $first_set->[$i][1] - $second_set->[$i][1] )**2 +
                 ( $first_set->[$i][2] - $second_set->[$i][2] )**2
    }

    # Devides by the number of member of the set.
    $rmsd = $rmsd / scalar @{ $first_set };

    return sqrt $rmsd;
}

#
# Calculates side-chain RMSD combinations.
# Input:
#     ${first,second}_atom_site - atom site data structure;
#     $unique_residue_key - unique residue key;
#     $options{'best_case'} - chooses those rmsd that are lowest.
# Output:
#     @sidechain_comparison_data - list of RMSD comparison data.
#

sub rmsd_sidechains
{
    my ( $parameters, $first_atom_site, $second_atom_site, $unique_residue_key,
         $options ) = @_;
    my ( $best_case ) = ( $options->{'best_case'} );

    $best_case //= 0;

    my $sig_figs_max = $parameters->{'_[local]_constants'}{'sig_figs_max'};
    # TODO: think if using of symmetric atom data should be optional or
    # mandatory.
    my $symmetrical_atom_names =$parameters->{'_[local]_symmetrical_atom_names'};

    my ( $residue_id, $chain, $pdbx_model_num ) = split /,/,$unique_residue_key;

    my @first_alt_ids =
        uniq @{ filter( { 'atom_site' => $first_atom_site,
                          'data' => [ 'label_alt_id' ],
                          'is_list' => 1 } ) };
    my @second_alt_ids =
        uniq @{ filter( { 'atom_site' => $second_atom_site,
                          'data' => [ 'label_alt_id' ],
                          'is_list' => 1 } ) };

    my @sidechain_comparison_data = ();
    for my $first_alt_id ( @first_alt_ids ) {
        my $first_sidechain_data =
            filter( { 'atom_site' => $first_atom_site,
                      'include' =>
                          { 'label_seq_id' => [ $residue_id ],
                            'label_asym_id' => [ $chain ],
                            'pdbx_PDB_model_num' => [ $pdbx_model_num ],
                            'label_alt_id' => [ $first_alt_id ] },
                      'data' =>
                          [ '[local]_selection_group', 'id',
                            'label_atom_id', 'label_seq_id',
                            'label_comp_id', 'label_asym_id',
                            'pdbx_PDB_model_num', 'label_alt_id',
                            'Cartn_x', 'Cartn_y', 'Cartn_z' ] } );
        $first_sidechain_data =
            [ sort { $a->[2] cmp $b->[2] } @{ $first_sidechain_data } ];

        my $residue_name = $first_sidechain_data->[0][4];
        my $rmsd_average;
        for my $second_alt_id ( @second_alt_ids ) {
            my $second_sidechain_data =
                filter( { 'atom_site' => $second_atom_site,
                          'include' =>
                              { 'label_seq_id' => [ $residue_id ],
                                'label_asym_id' => [ $chain ],
                                'pdbx_PDB_model_num' => [ $pdbx_model_num ],
                                'label_alt_id' => [ $second_alt_id ] },
                          'data' =>
                              [ 'id', '[local]_selection_group', 'label_atom_id',
                                'label_seq_id', 'label_comp_id', 'label_asym_id',
                                'pdbx_PDB_model_num', 'label_alt_id',
                                'Cartn_x', 'Cartn_y', 'Cartn_z' ] } );
            $second_sidechain_data =
                [ sort { $a->[2] cmp $b->[2] } @{ $second_sidechain_data } ];

            my %second_sidechain = ();
            for my $second_atom_data ( @{ $second_sidechain_data } ) {
                $second_sidechain{$second_atom_data->[2]} = $second_atom_data;
            }

            # Checks the length of the atom sets.
            # TODO: error message is duplicated in Measure::rmsd().
            if( scalar @{ $first_sidechain_data } ne
                scalar @{ $second_sidechain_data } ) {
                confess 'comparing different sizes of sets of the atoms ' .
                    'is not allowed';
            }

            my @current_sidechain_data = ();
            for( my $i = 0; $i <= $#{ $first_sidechain_data }; $i++ ) {
                if( $first_sidechain_data->[$i][2] ne
                    $second_sidechain_data->[$i][2] ) {
                    confess 'atom names do not match: ' .
                            "$first_sidechain_data->[$i][2] and " .
                            "$second_sidechain_data->[$i][2]";
                }

                # Checks if the side-chain have structural symmetry and chooses
                # the best RMSD value.
                my $symmetrical_atom_names =
                    $symmetrical_atom_names->{$residue_name}
                                             {$first_sidechain_data->[$i][2]};
                # TODO: make sure that there are no situations were best case
                # rmsd does not produce higher avarage rmsd.
                if( defined $symmetrical_atom_names ) {
                    my $rmsd;
                    my $symmetric_atom_data;
                    my @second_atom_names = ( $first_sidechain_data->[$i][2],
                                              @{ $symmetrical_atom_names } );
                    for my $second_atom_name ( @second_atom_names ) {
                        my $current_symmetric_atom_data =
                            $second_sidechain{$second_atom_name};
                        my $current_rmsd =
                            rmsd([[(@{$first_sidechain_data->[$i]})[8..10]]],
                                 [[(@{$current_symmetric_atom_data})[8..10]]]);
                        if( ! defined $rmsd || $current_rmsd < $rmsd ){
                            $symmetric_atom_data = $current_symmetric_atom_data;
                            $rmsd = $current_rmsd;
                        }
                    }

                    push @current_sidechain_data,
                        [ (@{$first_sidechain_data->[$i]})[0..7],
                          (@{$symmetric_atom_data})[0..7],
                          sprintf $sig_figs_max, $rmsd ];
                } else {
                    my $rmsd =
                        rmsd([[(@{$first_sidechain_data->[$i]})[8..10]]],
                             [[(@{$second_sidechain_data->[$i]})[8..10]]]);
                    push @current_sidechain_data,
                        [ (@{$first_sidechain_data->[$i]})[0..7],
                          (@{$second_sidechain_data->[$i]})[0..7],
                          sprintf $sig_figs_max, $rmsd ];
                }
            }

            if( $best_case ) {
                my $current_rmsd_average =
                    sum( map { $_->[-1] } @current_sidechain_data ) /
                    scalar @current_sidechain_data;
                if( ( ! @sidechain_comparison_data && ! defined $rmsd_average )||
                    $current_rmsd_average < $rmsd_average ) {
                    @sidechain_comparison_data = @current_sidechain_data;
                    $rmsd_average = $current_rmsd_average;
                }
            } else {
                push @sidechain_comparison_data, @current_sidechain_data;
            }
        }
    }

    return \@sidechain_comparison_data;
}

#
# Calculates energy values of the given structure.
# Input:
# Output:
#

sub energy
{
    my ( $parameters, $atom_site, $potential, $options  ) = @_;
    my ( $target_atom_ids, $only_sidechains, $decompose, $pairwise ) = (
        $options->{'target_atom_ids'},
        $options->{'only_sidechains'},
        $options->{'decompose'},
        $options->{'pairwise'},
    );

    $target_atom_ids //= [ sort keys %{ $atom_site } ];
    $only_sidechains //= 0;
    $decompose //= 0;

    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};
    my $interaction_atom_names =
        $parameters->{'_[local]_interaction_atom_names'};

    # Splits atom site into groups by its uniqueness.
    my $atom_site_groups = split_by( { 'atom_site' => $atom_site,
                                       'attributes' => [ 'pdbx_PDB_model_num',
                                                         'label_alt_id',
                                                         'label_asym_id' ],
                                       'append_dot' => 1 } );

    my $calculation_id = 1;
    my %energy = ();

    my %bonded_potentials = ();
    if( exists $potentials{$potential}{'bonded'} ) {
        %bonded_potentials = %{ $potentials{$potential}{'bonded'} };
    }

    my %non_bonded_potentials = ();
    if( exists $potentials{$potential}{'non_bonded'} ) {
        %non_bonded_potentials = %{ $potentials{$potential}{'non_bonded'} };
    }

    my %options = ();
    $options{'atom_site'} = $atom_site;

    my ( $grid_box, $target_cells ) = grid_box( $parameters,
                                                $atom_site,
                                                $edge_length_interaction,
                                                $target_atom_ids );

    my $neighbour_cells = identify_neighbour_cells( $grid_box, $target_cells );

    my @residue_energies = ();
    for my $cell ( sort { $a cmp $b } keys %{ $target_cells } ) {
        for my $atom_id ( @{ $target_cells->{$cell} } ) {
            next if ( any { $atom_site->{$atom_id}{'label_atom_id'} eq $_ }
                         @{ $interaction_atom_names } ) && $only_sidechains;

            my @residue_energy = ();

            # Calculates bonded potential energy term.
            my %bonded_residue_energy = (); # For faster neighbour energy search.
            for my $bonded_potential ( keys %bonded_potentials ) {
                my $residue_bonded_energy =
                    $bonded_potentials{$bonded_potential}( $parameters,
                                                           $atom_id,
                                                           \%options );
                for my $bonded_energy ( @{ $residue_bonded_energy } ) {
                    my $neighbour_atom_id = $bonded_energy->atoms->[-1];
                    push @{ $bonded_residue_energy{$atom_id}{$neighbour_atom_id} },
                        $bonded_energy;
                }
            }

            for my $neighbour_atom_id ( uniq @{ $neighbour_cells->{$cell} } ) {
                if( ( $atom_id ne $neighbour_atom_id ) &&
                    ( ! is_neighbour( $atom_site,
                                      $atom_id,
                                      $neighbour_atom_id ) ) &&
                    ( ! is_second_neighbour( $atom_site,
                                             $atom_id,
                                             $neighbour_atom_id ) ) ) {

                    # Adds bonded potential energy term.
                    for my $bonded_potential (
                        @{ $bonded_residue_energy{$atom_id}{$neighbour_atom_id} } ) {
                        push @residue_energy, $bonded_potential;
                    }

                    # Adds non-bonded potential energy term.
                    for my $non_bonded_potential ( keys %non_bonded_potentials ) {
                        my $energy_potential = Energy->new();
                        $energy_potential->set_energy(
                            $non_bonded_potential,
                            [ $atom_id, $neighbour_atom_id ],
                            $non_bonded_potentials{$non_bonded_potential}(
                                $parameters,
                                $atom_site->{$atom_id},
                                $atom_site->{$neighbour_atom_id},
                                \%options
                            )
                        );

                        push @residue_energy, $energy_potential;
                    }
                }
            }

            push @residue_energies, @residue_energy;
        }
    }

    my @atom_pairs = ();
    my %atom_pair_interactions = ();
    for my $residue_energy ( @residue_energies ) {
        my $atom_id = $residue_energy->atoms->[0];
        my $interaction_atom_id = $residue_energy->atoms->[-1];
        my $residue_potential = $residue_energy->energy_type;

        if( ! exists $atom_pair_interactions{$atom_id}
                                            {$interaction_atom_id}
                                            {$residue_potential} ) {
            push @atom_pairs, [ $atom_id, $interaction_atom_id ];
        }

        push @{ $atom_pair_interactions{$atom_id}
                                       {$interaction_atom_id}
                                       {$residue_potential} },
            $residue_energy;
    }

    my @energies = ();
    for my $atom_pair ( @atom_pairs ) {
        my $atom_id = $atom_pair->[0];
        my $interaction_atom_id = $atom_pair->[1];
        my @energy_types =
            sort keys %{$atom_pair_interactions{$atom_id}{$interaction_atom_id}};

        if( $pairwise && $decompose ) {
            for my $energy_type ( @energy_types ) {
                push @energies,
                    @{ $atom_pair_interactions{$atom_id}
                                              {$interaction_atom_id}
                                              {$energy_type} };
            }
        } elsif( $decompose ) {

        } elsif( $pairwise ) {

        } else {
            my $energy_sum_value = 0;
            for my $energy_type ( @energy_types ) {
                $energy_sum_value +=
                    $atom_pair_interactions{$atom_id}
                                           {$interaction_atom_id}
                                           {$energy_type}[0]->value;
            }

            my $energy_sum = Energy->new();
            $energy_sum->set_energy( $potential,
                                     [ $atom_id, $interaction_atom_id ],
                                     $energy_sum_value );
            push @energies, $energy_sum;
        }
    }

    return \@energies;
}

1;
