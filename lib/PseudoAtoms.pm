package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( add_hydrogens
                     calc_favourable_angle
                     calc_favourable_angles
                     generate_library
                     generate_pseudo
                     generate_rotamer
                     library_to_csv
                     replace_with_rotamer );

use B qw( svref_2object );
use Carp;
use Clone qw( clone );
use List::Util qw( max );
use List::MoreUtils qw( any
                        uniq );
use Math::Trig qw( acos );
use threads;

use AtomInteractions qw( hard_sphere
                         soft_sphere
                         leonard_jones
                         composite );
use Combinatorics qw( permutation );
use ConnectAtoms qw( append_connections
                     connect_atoms
                     is_neighbour
                     is_second_neighbour );
use Constants qw( $EDGE_LENGTH_INTERACTION
                  $PI
                  $SIG_FIGS_MIN );
use ForceField::General;
use Grid qw( grid_box
             identify_neighbour_cells );
use LinearAlgebra qw( mult_matrix_product
                      switch_ref_frame );
use Measure qw( all_dihedral
                dihedral_angle
                bond_angle );
use BondProperties qw( hybridization
                       rotatable_bonds
                       unique_rotatables );
use Multithreading qw( multithreading );
use PDBxParser qw( create_pdbx_entry
                   determine_residue_keys
                   filter
                   filter_by_unique_residue_key
                   split_by
                   unique_residue_key );
use Sampling qw( sample_angles );
use SidechainModels qw( rotation_only );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates pseudo-atoms from side chain models that are written in list of
# analytical matrices.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser). Must have
#     any sidechain model function applied to it (see SidechainModels.pm);
#     $args->{atom_specifier} - hash of hashes for selecting atoms by attributes
#     (see PDBxParser.pm);
#     $args->{angle_values} - hash of arrays that describe possible values of
#     dihedral angles:
#     Ex.: { 'chi1' => [ 0, 0.4, 1.5, 2.0 ],
#            'chi2' => [ 0, 2 ] };
#     $args->{last_atom_id} - last atom id for assigning new ids for pseud
#     atoms;
#     $args->{alt_group_id} - alternative group id that is used to distinguish
#     pseudo atoms. Very useful when generating rotamers.
# Output:
#     $pseudo_atom_site - atom site data structure for pseudo-atoms with
#     additional 'conformation' attribute.
#

sub generate_pseudo
{
    my ( $args ) = @_;
    my ( $atom_site, $atom_specifier, $angle_values, $last_atom_id,
         $alt_group_id ) =
        ( $args->{'atom_site'}, $args->{'atom_specifier'},
          $args->{'angle_values'}, $args->{'last_atom_id'},
          $args->{'alt_group_id'} );

    $last_atom_id //= max( keys %{ $atom_site } );
    $alt_group_id //= 1;

    my %atom_site = %{ clone( $atom_site ) };
    my %pseudo_atom_site;

    my @atom_ids =
        @{ filter( { 'atom_site' => \%atom_site,
                     'include' => $atom_specifier,
                     'data' => [ 'id' ],
                     'is_list' => 1 } ) };

    for my $atom_id ( @atom_ids ) {
        my $conformation = $atom_site{"$atom_id"}{'conformation'};

        confess "atom with id $atom_id lacks 'conformation' key."
            if ! defined $conformation;

        # Calculates current dihedral angles of rotatable bonds. Will be used
        # for reseting dihedral angles to 0 degree angle.
        my $residue_unique_key = unique_residue_key( $atom_site{$atom_id} );

        my %angles =
            %{ all_dihedral(
                   filter_by_unique_residue_key( $atom_site,
                                                 $residue_unique_key, 1 ) ) };

        # Iterates through combinations of angles and evaluates conformational
        # model.
        my @angle_names = sort { $a cmp $b } keys %{ $angle_values };
        my @angle_values;
        for my $angle_name ( @angle_names ) {
            push @angle_values,
                 [ map { $_ - $angles{$residue_unique_key}
                                     {"$angle_name"}{'value'} }
                   @{ $angle_values->{"$angle_name"} } ];
        }

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
                { 'atom_site' => \%pseudo_atom_site,
                  'id' => $last_atom_id,
                  'type_symbol' => $atom_site{$atom_id}{'type_symbol'},
                  'label_atom_id' => $atom_site{$atom_id}{'label_atom_id'},
                  'label_alt_id' => $alt_group_id,
                  'label_comp_id' => $atom_site{$atom_id}{'label_comp_id'},
                  'label_asym_id' => $atom_site{$atom_id}{'label_asym_id'},
                  'label_entity_id' => $atom_site{$atom_id}{'label_entity_id'},
                  'label_seq_id' => $atom_site{$atom_id}{'label_seq_id'},
                  'cartn_x' => sprintf( $SIG_FIGS_MIN, $transf_atom_coord->[0][0] ),
                  'cartn_y' => sprintf( $SIG_FIGS_MIN, $transf_atom_coord->[1][0] ),
                  'cartn_z' => sprintf( $SIG_FIGS_MIN, $transf_atom_coord->[2][0] ),
                  'pdbx_PDB_model_num' =>
                      $atom_site{$atom_id}{'pdbx_PDB_model_num'},
                } );
            # Adds atom id that pseudo atoms was made of.
            $pseudo_atom_site{$last_atom_id}{'origin_atom_id'} = $atom_id;
            # Adds hybridization, connection, conformation data from origin atom.
            $pseudo_atom_site{$last_atom_id}{'hybridization'} =
                $atom_site{$atom_id}{'hybridization'};
            # FIXME: be careful - it might produce contradictions between new
            # and old pseudo atoms.
            $pseudo_atom_site{$last_atom_id}{'connections'} =
                $atom_site{$atom_id}{'connections'};
            $pseudo_atom_site{$last_atom_id}{'conformation'} =
                $atom_site{$atom_id}{'conformation'};
            # Adds information about used dihedral angle values and names.
            $pseudo_atom_site{$last_atom_id}{'dihedral_names'} = \@angle_names;
            $pseudo_atom_site{$last_atom_id}{'dihedral_angles'} =
                { map { ( $_ => $angle_values{$_} +
                                $angles{$residue_unique_key}{$_}{'value'} ) }
                  @angle_names };
            # Adds additional pseudo-atom flag for future filtering.
            $pseudo_atom_site{$last_atom_id}{'is_pseudo_atom'} = 1;
        }
    }

    return \%pseudo_atom_site;
}

#
# Generates rotamers according to given angle values.
# Input:
#     $args->{'atom_site'} - atom site data structure (see PDBxParser);
#     $args->{'angle_values'} - name and value of angles in hash form.
#     $args->{'last_atom_id'} - last atom id for assigning new ids for pseudo atoms;
#     $args->{'alt_group_id'} - alternative group id that is used to distinguish
#     pseudo atoms;
#     $args->{'set_missing_angles_to_zero'} - if angles are unknown, they are
#     set to 0 rad;
#     $args->{'keep_origin_id'} - keeps ids the same as original atom's. Very
#     useful when replacing existing residues.
# Output:
#     %generated_rotamers - atom site data structure with additional rotamer
#     data.
#

sub generate_rotamer
{
    my ( $args ) = @_;
    my ( $atom_site, $angle_values, $last_atom_id, $alt_group_id,
         $set_missing_angles_to_zero, $keep_origin_id, $keep_origin_alt_id ) =
        ( $args->{'atom_site'}, $args->{'angle_values'}, $args->{'last_atom_id'},
          $args->{'last_atom_id'}, $args->{'set_missing_angles_to_zero'},
          $args->{'keep_origin_id'}, $args->{'keep_origin_alt_id'} );

    $last_atom_id //= max( keys %{ $atom_site } );
    $alt_group_id //= 1;
    $set_missing_angles_to_zero //= 0;
    $keep_origin_id //= 0;
    $keep_origin_alt_id //= 0; # Has higher priority than $alt_group_id.

    my %atom_site = %{ clone( $atom_site ) };
    my %rotamer_atom_site;

    for my $residue_unique_key ( keys %{ $angle_values } ) {
        my ( $residue_id, $residue_chain, $pdbx_model, $residue_alt_id ) =
            split /,/, $residue_unique_key;
        $residue_alt_id =
            $keep_origin_alt_id ? $residue_alt_id : $alt_group_id;
        my $residue_site =
            filter_by_unique_residue_key( \%atom_site, $residue_unique_key, 1 );

        my $rotatable_bonds = rotatable_bonds( $residue_site );

        for my $atom_id ( sort { $a <=> $b } keys %{ $residue_site } ) {
            if( ! exists $rotatable_bonds->{$atom_id} ) { next; }

            my %angles;
            for my $angle_name ( keys %{ $rotatable_bonds->{$atom_id} } ) {
                if( exists $angle_values->{"$residue_unique_key"}{$angle_name} &&
                    defined $angle_values->{"$residue_unique_key"}{$angle_name}){
                    $angles{$angle_name} =
                        [ $angle_values->{"$residue_unique_key"}{$angle_name} ];
                } else {
                    if( $set_missing_angles_to_zero ) {
                        $angles{$angle_name} = [ 0.0 ];
                    } else {
                        confess "no values for $angle_name were assigned.";
                    }
                }
            }

            %rotamer_atom_site =
                ( %rotamer_atom_site,
                  %{ generate_pseudo( {
                         'atom_site' => { ( %atom_site, %rotamer_atom_site ) },
                         'atom_specifier' => { 'id' => [ $atom_id ] },
                         'angle_values' => \%angles,
                         'last_atom_id' => $last_atom_id,
                         'alt_group_id' => $residue_alt_id } ) } );
            $last_atom_id++;
        }
    }

    if( $keep_origin_id ) {
        my %rotamer_atom_site_old_ids;
        for my $atom_id ( keys %rotamer_atom_site ) {
            my $origin_id = $rotamer_atom_site{$atom_id}{'origin_atom_id'};
            $rotamer_atom_site{$atom_id}{'id'} = $origin_id;
            $rotamer_atom_site_old_ids{$origin_id}= $rotamer_atom_site{$atom_id};
        }
        return \%rotamer_atom_site_old_ids;
    }

    return \%rotamer_atom_site;
}

#
# Generates rotamer libraries by specified arguments that include atom
# movements and interactions between atoms.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $args->{residue_unique_keys} - array of unique residue keys
#     (see PDBxParser::unique_residue_key);
#     $args->{include_interactions} - selection data structure
#     (see PDBxParser::filter) that is used to select atoms that will be included
#     into calculations of energy;
#     $args->{small_angle} - angle by which rotation is made;
#     $args->{conf_model} - possible sidechain movements described by sidechain
#     model functions in SidechainModels.pm;
#     $args->{interactions} - interaction models described by functions in
#     AtomInteractions.pm;
#     $args->{parameters} - parameters that are passed to interaction function;
#     $args->{energy_cutoff_atom} - maximum amount of energy that is allowed for
#     atom to have in the rotamer according to potential function;
#     $args->{energy_cutoff_residue} - maximum amount of energy (the sum of
#     energies of all atoms) that is allowed for residue to have in the rotamer
#     according to potential function;
#     $args->{threads} - number of threads.
# Output:
#     %library_atom_site - atom site data structure with additional data.
#

sub generate_library
{
    my ( $args ) = @_;
    my $atom_site = $args->{'atom_site'};
    my $residue_unique_keys = $args->{'residue_unique_keys'};
    my $include_interactions = $args->{'include_interactions'};
    my $small_angle = $args->{'small_angle'};
    my $conf_model = $args->{'conf_model'};
    my $interactions = $args->{'interactions'};
    my $parameters = $args->{'parameters'};
    my $energy_cutoff_atom = $args->{'energy_cutoff_atom'};
    my $energy_cutoff_residue = $args->{'energy_cutoff_residue'};
    my $is_hydrogen_explicit = $args->{'is_hydrogen_explicit'};
    my $threads = $args->{'threads'};

    $conf_model //= 'rotation_only';
    $energy_cutoff_residue //= 'Inf';
    $threads //= 1;
    $include_interactions //= { 'label_atom_id' =>
                                    \@General::INTERACTION_ATOM_NAMES };
    $is_hydrogen_explicit //= 0;

    # Selection of potential function.
    my %potential_functions = ( 'composite' => \&composite,
                                'hard_sphere' => \&hard_sphere,
                                'soft_sphere' => \&soft_sphere,
                                'leonard_jones' => \&leonard_jones, );
    my $potential_function = $potential_functions{"$interactions"};

    my %rotamer_library;

    my $atom_site_groups =
        split_by( { 'atom_site' => $atom_site,
                    'attributes' => [ 'pdbx_PDB_model_num', 'label_alt_id' ],
                    'append_dot' => 1  } );

    for my $atom_site_identifier ( sort keys %{ $atom_site_groups } ) {
        my ( $pdbx_model_num, $alt_id ) = split /,/, $atom_site_identifier;
        my $current_atom_site =
            filter( { 'atom_site' => $atom_site,
                      'include' =>
                          {'id' => $atom_site_groups->{$atom_site_identifier}}});

        # Predetermines geometrically clearly defined hydrogen positions.
        my $current_atom_site_w_H = clone( $current_atom_site );
        my $hydrogens =
            add_hydrogens( $current_atom_site_w_H,
                           { 'alt_group_id' => q{.},
                             'add_only_clear_positions' => 1,
                             'use_origins_alt_group_id' => 1 } );
        append_connections( $current_atom_site_w_H, $hydrogens );
        $current_atom_site_w_H = { %{ $current_atom_site_w_H }, %{ $hydrogens }};
        hybridization( $current_atom_site_w_H );

        # Also, prepares hydrogen-free structure.
        my $current_atom_site_no_H = clone(
            filter( { 'atom_site' => $current_atom_site_w_H,
                      'exclude' => { 'type_symbol' => [ 'H' ] } } )
        );
        connect_atoms( $current_atom_site_no_H );
        hybridization( $current_atom_site_no_H );

        # Generates conformational models before checking for
        # clashes/interactions for given residues.
        if( $conf_model eq 'rotation_only' ) {
            rotation_only( $current_atom_site_no_H );
            rotation_only( $current_atom_site_w_H );
        } else {
            confess 'conformational model was not defined.';
        }

        # Finds where CA of target residues are.
        my @target_ca_ids;
        for my $residue_unique_key ( @{ $residue_unique_keys } ) {
            my $residue_site =
                filter_by_unique_residue_key( $current_atom_site,
                                              $residue_unique_key, 1 );
            my $atom_ca_id =
                filter( { 'atom_site' => $residue_site,
                          'include' => { 'label_atom_id' => [ 'CA' ] },
                          'data' => [ 'id' ],
                          'is_list' => 1 } )->[0];
            push @target_ca_ids, $atom_ca_id;
        }

        # Creates the grid box that has edge length of sum of all bonds of the
        # longest side-chain branch in arginine. Length: 3 * (C-C) + (C-N) + 2
        # * (C=N) + (N-H).
        # TODO: should consider using shorter distances, because bonds have limits
        # on maximum bending and having shorter edge length reduces calculation
        # time.
        my ( $grid_box, $target_cell_idxs ) =
            grid_box( $current_atom_site_w_H, $EDGE_LENGTH_INTERACTION,
                      \@target_ca_ids );

        my $neighbour_cells =
            identify_neighbour_cells( $grid_box, $target_cell_idxs );

        for my $cell ( sort { $a cmp $b } keys %{ $target_cell_idxs } ) {
            for my $ca_atom_id ( @{ $target_cell_idxs->{$cell} } ) {
                my $residue_id =
                    $current_atom_site->{$ca_atom_id}{'label_seq_id'};
                my $residue_chain =
                    $current_atom_site->{$ca_atom_id}{'label_asym_id'};
                my $residue_site =
                    filter( { 'atom_site' => $current_atom_site,
                              'include' =>
                                  { 'pdbx_PDB_model_num' => [ $pdbx_model_num ],
                                    'label_alt_id' => [ $alt_id, '.' ],
                                    'label_seq_id' => [ $residue_id ],
                                    'label_asym_id' => [ $residue_chain ] } } );
                my $residue_unique_key =
                    determine_residue_keys( $residue_site,
                                            {'exclude_dot' => 1} )->[0];

                # Because the change of side-chain position might impact the
                # surrounding, iteraction site consists of only main chain atoms.
                my %interaction_site_no_H =
                    %{ filter( { 'atom_site' => $current_atom_site_no_H,
                                 'include' =>
                                     { 'id' => $neighbour_cells->{$cell},
                                       %{ $include_interactions } } } ) };

                # First, checks angles by step-by-step adding atoms to sidechains.
                # This is called growing side chain.
                my @allowed_angles =
                    @{ calc_favourable_angles(
                           { 'atom_site' => $current_atom_site_no_H,
                             'residue_unique_key' => $residue_unique_key,
                             'interaction_site' => \%interaction_site_no_H,
                             'small_angle' => $small_angle,
                             'potential_function' => $potential_function,
                             'energy_cutoff_atom' => $energy_cutoff_atom,
                             'parameters' => $parameters,
                             'threads' => $threads } ) };

                # Then, re-checks if each atom of the rotamer obey energy
                # cutoffs.
                my %interaction_site_w_H =
                    %{ filter( { 'atom_site' => $current_atom_site_w_H,
                                 'include' =>
                                     { 'id' => $neighbour_cells->{$cell},
                                       %{ $include_interactions } } } ) };

                my ( $allowed_angles, $energy_sums ) =
                    @{ calc_full_atom_energy(
                           { 'atom_site' => ( $is_hydrogen_explicit ?
                                              $current_atom_site_w_H :
                                              $current_atom_site_no_H ),
                             'residue_unique_key' => $residue_unique_key,
                             'interaction_site' => ( $is_hydrogen_explicit ?
                                                     \%interaction_site_w_H :
                                                     \%interaction_site_no_H ),
                             'small_angle' => $small_angle,
                             'potential_function' => $potential_function,
                             'energy_cutoff_atom' => $energy_cutoff_atom,
                             'is_hydrogen_explicit' => $is_hydrogen_explicit,
                             'parameters' => $parameters },
                           [ @allowed_angles ] ) };

                # my ( $allowed_angles, $energy_sums ) =
                #     @{ multithreading(
                #            \&calc_full_atom_energy,
                #            { 'atom_site' => $current_atom_site_w_H,
                #              'residue_unique_key' => $residue_unique_key,
                #              'interaction_site' => \%interaction_site_w_H,
                #              'small_angle' => $small_angle,
                #              'potential_function' => $potential_function,
                #              'energy_cutoff_atom' => $energy_cutoff_atom,
                #              'is_hydrogen_explicit' => $is_hydrogen_explicit,
                #              'parameters' => $parameters },
                #            [ @allowed_angles ],
                #            $threads ) };

                if( ! @{ $allowed_angles } ) {
                    confess "no possible rotamer solutions were detected.\n";
                }

                for( my $i = 0; $i <= $#{ $allowed_angles }; $i++  ) {
                    my %angles =
                        map { my $angle_id = $_ + 1;
                              ( "chi$angle_id" => $allowed_angles->[$i][$_])}
                            ( 0..$#{ $allowed_angles->[$i] } );
                    my $rotamer_energy_sum = $energy_sums->[$i];
                    if( defined $rotamer_energy_sum &&
                        $rotamer_energy_sum <= $energy_cutoff_residue ) {
                        push @{ $rotamer_library{"$residue_unique_key"} },
                            { 'angles' => \%angles,
                              'potential' => $interactions,
                              'potential_energy_value' => $rotamer_energy_sum };
                    }
                }
            }
        }
    }

    return \%rotamer_library;
}

#
# Calculates favourable rotamer angles for a given residue.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $args->{residue_unique_key} - unique residue key
#     (see PDBxParser::unique_residue_key);
#     $args->{interaction_site} - atom data structure that is included into
#     energy calculations;
#     $args->{small_angle} - angle by which rotation is made;
#     $args->{potential_function} - reference to the potential function that is
#     used for calculating energy;
#     $args->{energy_cutoff_atom} - maximum amount of energy that is allowed for
#     atom to have in the rotamer according to potential function;
#     $args->{parameters} - parameters that are passed to interaction function;
#     $args->{threads} - number of threads.
# Output:
#     @allowed_angles - list of groups of allowed angles.
#     Ex.: ( [ 0.00, 3.14 ], [ 3.14, 6.28 ] ).
#

sub calc_favourable_angles
{
    my ( $args ) = @_;

    my ( $atom_site, $residue_unique_key, $interaction_site, $small_angle,
         $potential_function, $energy_cutoff_atom, $parameters, $threads ) = (
        $args->{'atom_site'},
        $args->{'residue_unique_key'},
        $args->{'interaction_site'},
        $args->{'small_angle'},
        $args->{'potential_function'},
        $args->{'energy_cutoff_atom'},
        $args->{'parameters'},
        $args->{'threads'},
    );

    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $residue_unique_key, 1 );

    my $rotatable_bonds = rotatable_bonds( $residue_site );
    if( ! %{ $rotatable_bonds } ) { next; }

    # Goes through each atom in side chain and calculates interaction
    # potential with surrounding atoms. CA and CB are non-movable atoms
    # so, they are marked as starting atoms.
    my $ca_atom_id =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'CA' ] },
                  'data' => [ 'id' ] } )->[0][0];
    my $cb_atom_id =
        filter( { 'atom_site' => $residue_site,
                  'include' => { 'label_atom_id' => [ 'CB' ] },
                  'data' => [ 'id' ] } )->[0][0];

    my @visited_atom_ids = ( $ca_atom_id, $cb_atom_id );
    my @next_atom_ids =
        grep { $_ ne $ca_atom_id }
            @{ $residue_site->{$cb_atom_id}{'connections'} };

    my @sampled_angles =
        map { [ $_ ] } @{ sample_angles( [ [ 0, 2 * $PI ] ], $small_angle ) };
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
            if( scalar( @{ $allowed_angles[0] } ) <
                scalar( keys %{ $rotatable_bonds->{$atom_id} } ) ) {
                @allowed_angles =
                    @{ permutation( 2, [], [ \@allowed_angles,
                                             \@sampled_angles ], [] ) };
                @allowed_energies =
                    @{ permutation( 2, [], [ \@allowed_energies,
                                             \@zero_energies ], [] ) };
                # Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
                @allowed_angles =
                    map { [ @{ $_->[0] }, @{ $_->[1] } ] } @allowed_angles;
                @allowed_energies =
                    map { [ $_->[0][0] ] } @allowed_energies;
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            # Marks neighbouring atoms.
            push @neighbour_atom_ids,
                 @{ $atom_site->{$atom_id}{'connections'} };

            # Starts calculating potential energy.
            my ( $next_allowed_angles, $next_allowed_energies ) =
                @{ calc_favourable_angle(
                       { 'atom_site' => $atom_site,
                         'atom_id' => $atom_id,
                         'interaction_site' => $interaction_site,
                         'energy_cutoff_atom' => $energy_cutoff_atom,
                         'potential_function' => $potential_function,
                         'parameters' => $parameters },
                       [ \@allowed_angles, \@allowed_energies, ] ) };

            # my ( $next_allowed_angles, $next_allowed_energies ) =
            #     @{ multithreading(
            #            \&calc_favourable_angle,
            #            { 'atom_site' => $atom_site,
            #              'atom_id' => $atom_id,
            #              'interaction_site' => $interaction_site,
            #              'energy_cutoff_atom' => $energy_cutoff_atom,
            #              'potential_function' => $potential_function,
            #              'parameters' => $parameters },
            #            [ \@allowed_angles, \@allowed_energies, ],
            #            $threads ) };

            if( scalar @{ $next_allowed_angles } > 0 ) {
                @allowed_angles = @{ $next_allowed_angles };
                @allowed_energies = @{ $next_allowed_energies };
            } else {
                confess "no possible rotamer solutions were detected.\n";
            }
        }

        # Determines next atoms that should be visited.
        @next_atom_ids = (); # Resets value for the new ones to be
                             # appended.
        for my $neighbour_atom_id ( uniq @neighbour_atom_ids ) {
            if( ( ! any { $neighbour_atom_id eq $_ } @visited_atom_ids ) ) {
                push @next_atom_ids, $neighbour_atom_id;
            }
        }
    }

    return \@allowed_angles;
}

#
# Calculates energy values for given rotamer angles.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $args->{atom_id} - atom id;
#     $args->{interaction_site} - atom data structure that is included into
#     energy calculations;
#     $args->{potential_function} - reference to the potential function that is
#     used for calculating energy;
#     $args->{energy_cutoff_atom} - maximum amount of energy that is allowed for
#     atom to have in the rotamer according to potential function;
#     $args->{parameters} - parameters that are passed to interaction function.
# Output:
#     @allowed_angles - list of groups of allowed angles.
#     Ex.: ( [ 0.00, 3.14 ], [ 3.14, 6.28 ] );
#     @allowed_energies - list of energies of the allowed angles.
#     Ex.: ( [ -2.32 ], [ -15.01 ] ).
#

sub calc_favourable_angle
{
    my ( $args, $array_blocks ) = @_;

    my ( $atom_site, $atom_id, $interaction_site, $potential_function,
         $energy_cutoff_atom, $parameters ) = (
        $args->{'atom_site'},
        $args->{'atom_id'},
        $args->{'interaction_site'},
        $args->{'potential_function'},
        $args->{'energy_cutoff_atom'},
        $args->{'parameters'},
    );

    my %parameters = defined $parameters ? %{ $parameters } : ();
    $parameters{'atom_site'} = $atom_site;

    my @allowed_angles;
    my @allowed_energies;
    for( my $i = 0; $i <= $#{ $array_blocks->[0] }; $i++ ) {
        my $angles = $array_blocks->[0][$i];
        my $energies = $array_blocks->[1][$i][0];
        my %angles =
            map { my $angle_id = $_ + 1; ( "chi$angle_id" => [ $angles->[$_] ] )}
                0..$#{ $angles };

        my $pseudo_atom_site =
            generate_pseudo( { 'atom_site' => $atom_site,
                               'atom_specifier' => { 'id' => [ "$atom_id" ] },
                               'angle_values' => \%angles } );
        my $pseudo_atom_id = ( keys %{ $pseudo_atom_site } )[0];
        my $pseudo_origin_id =
            $pseudo_atom_site->{$pseudo_atom_id}{'origin_atom_id'};

        my $potential_energy = 0; # TODO: look if here should be zeros.
        my $potential_sum = 0;
        foreach my $interaction_id ( keys %{ $interaction_site } ) {
            if( ( ! is_neighbour( $atom_site,
                                  $pseudo_origin_id,
                                  $interaction_id ) ) &&
                ( ! is_second_neighbour( $atom_site,
                                         $pseudo_origin_id,
                                         $interaction_id ) ) ) {
                $potential_energy =
                    $potential_function->(
                        $pseudo_atom_site->{$pseudo_atom_id},
                        $atom_site->{$interaction_id},
                        \%parameters );
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
            push @allowed_angles, $angles;
            push @allowed_energies, [ $energies + $potential_sum ];
        }
    }

    return [ \@allowed_angles, \@allowed_energies ];
}

#
# Calculates full atom energy by including even the atoms of the current rotamer.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $args->{residue_unique_keys} - array of unique residue keys
#     (see PDBxParser::unique_residue_key);
#     $args->{interaction_site} - atom data structure that is included into
#     energy calculations;
#     $args->{small_angle} - angle by which rotation is made;
#     $args->{potential_function} - reference to the potential function that is
#     used for calculating energy;
#     $args->{energy_cutoff_atom} - maximum amount of energy that is allowed for
#     atom to have in the rotamer according to potential function;
#     $args->{parameters} - parameters that are passed to interaction function;
#     $array_blocks - an array of arrays that contain suggested rotamer angles.
# Output:
#     @allowed_angles - list of groups of allowed angles.
#     Ex.: ( [ 0.00, 3.14 ], [ 3.14, 6.28 ] );
#     @energy_sums - list of energy sums of the all atoms in the residue after
#     bond is rotated according to @allowed_angles.
#     Ex.: ( [ -45.02 ], [ -15.00 ] ).
#

sub calc_full_atom_energy
{
    my ( $args, $array_blocks ) = @_;

    my ( $atom_site, $residue_unique_key, $interaction_site, $small_angle,
         $potential_function, $energy_cutoff_atom, $is_hydrogen_explicit,
         $parameters ) = (
        $args->{'atom_site'},
        $args->{'residue_unique_key'},
        $args->{'interaction_site'},
        $args->{'small_angle'},
        $args->{'potential_function'},
        $args->{'energy_cutoff_atom'},
        $args->{'is_hydrogen_explicit'},
        $args->{'parameters'},
    );

    # Creates all-atom model (even with uncertain hydrogen positions) for
    # selected residue.
    my ( $residue_id, $residue_chain, $pdbx_model_num,
         $residue_alt ) = split /,/sxm, $residue_unique_key;
    my $residue_site =
        filter( { 'atom_site' => $atom_site,
                  'include' => { 'label_seq_id' => [ $residue_id ],
                                 'label_asym_id' => [ $residue_chain ],
                                 'label_alt_id' => [ $residue_alt, '.' ],
                                 'pdbx_PDB_model_num' => [ $pdbx_model_num ] }});

    # Adds hydrogens to the residue_site.
    if( $is_hydrogen_explicit ) {
        my $hydrogens =
            add_hydrogens( $residue_site,
                           { 'use_existing_connections' => 1,
                                 'use_existing_hybridizations' => 1,
                                 'exclude_by_atom_name' => [ 'N', 'C' ],
                                 'use_origins_alt_group_id' => 1 } );
        append_connections( $residue_site, $hydrogens );
        $residue_site = { %{ $residue_site }, %{ $hydrogens } };

        rotation_only( $residue_site );
    }

    # Identifies missing unique rotatable bonds.
    my %uniq_rotatable_bonds = %{ unique_rotatables( $residue_site ) };

    if( ! %uniq_rotatable_bonds ) { next; }

    my $uniq_rotatable_bond_num = scalar keys %uniq_rotatable_bonds;
    my $missing_rotatable_bond_num = $uniq_rotatable_bond_num -
                                     scalar @{ $array_blocks->[0] };

    my @checkable_angles = @{ $array_blocks };
    if( $missing_rotatable_bond_num ) {
        my @sampled_angles =
            map { [ $_ ] } @{sample_angles( [ [ 0, 2 * $PI ] ], $small_angle )};
        foreach( 1..$missing_rotatable_bond_num ) {
            @checkable_angles =
                @{ permutation( 2, [], [ \@checkable_angles,
                                         \@sampled_angles ], [] ) };
            # Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
            @checkable_angles =
                map { [ @{ $_->[0] }, @{ $_->[1] } ] } @checkable_angles;
        }
    }

    # Checks for inter-atom interactions and determines if energies
    # comply with cutoffs.
    my @allowed_angles;
    my @energy_sums;
  ALLOWED_ANGLES:
    for( my $i = 0; $i <= $#checkable_angles; $i++ ) {
        my %angles =
            map { my $angle_id = $_ + 1; ( "chi$angle_id" =>
                                               $checkable_angles[$i][$_] ) }
                0..$#{ $checkable_angles[$i] };

        my %rotamer_site = %{ $residue_site };
        replace_with_rotamer( \%rotamer_site, $residue_unique_key, \%angles );

        my @rotamer_atom_ids =
            sort keys %{ filter( { 'atom_site' => \%rotamer_site,
                                   'exclude' =>
                                   { 'label_atom_id' =>
                                         \@General::INTERACTION_ATOM_NAMES } } ) };
        # HACK: make sure that $interaction_site atom ids are updated by
        # %rotamer_site
        my %rotamer_interaction_site = ( %{ $interaction_site }, %rotamer_site );

        # HACK: should connect_atoms() be used here?
        # connect_atoms( \%rotamer_interaction_site );

        $parameters->{'atom_site'} = \%rotamer_interaction_site;

        my $rotamer_energy_sum = 0;
        for my $rotamer_atom_id ( @rotamer_atom_ids ) {
            for my $neighbour_atom_id ( sort keys %rotamer_interaction_site ) {
                my $rotamer_atom_energy = 0;
                if( ( $rotamer_atom_id ne $neighbour_atom_id ) &&
                    ( ! is_neighbour( \%rotamer_interaction_site,
                                      $rotamer_atom_id,
                                      $neighbour_atom_id ) ) &&
                    ( ! is_second_neighbour( \%rotamer_interaction_site,
                                             $rotamer_atom_id,
                                             $neighbour_atom_id ) ) ){
                    $rotamer_atom_energy +=
                        $potential_function->(
                            $rotamer_interaction_site{$rotamer_atom_id},
                            $rotamer_interaction_site{$neighbour_atom_id},
                            $parameters );

                    next ALLOWED_ANGLES
                        if $rotamer_atom_energy > $energy_cutoff_atom;

                    $rotamer_energy_sum += $rotamer_atom_energy;
                }
            }
        }

        push @allowed_angles, $checkable_angles[$i];
        push @energy_sums, $rotamer_energy_sum;
    }

    return [ \@allowed_angles, \@energy_sums ] ;
}

#
# Rotates residue bonds by specified dihedral angles.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $residue_unique_key - unique residue key
#     (see PDBxParser::unique_residue_key);
#     $angle_values - name and value of angles in hash form.
# Output:
#     changes coordinates of selected residue due to bond rotations.
#

sub replace_with_rotamer
{
    my ( $atom_site, $residue_unique_key, $angle_values ) = @_;

    my ( undef, undef, undef, $alt_group_id ) = split /,/, $residue_unique_key;
    my $residue_site =
        generate_rotamer( { 'atom_site' => $atom_site,
                            'angle_values' =>
                                { $residue_unique_key => $angle_values  },
                            'alt_group_id' => 'X', # HACK: $keep_origin_alt_id
                                                   # should be used.
                            'set_missing_angles_to_zero' => 1 } );

    for my $residue_atom_id ( keys %{ $residue_site } ) {
        my $residue_origin_atom_id = $residue_site->{$residue_atom_id}
                                                    {'origin_atom_id'};
        $residue_site->{$residue_atom_id}{'id'} = $residue_origin_atom_id;
        $residue_site->{$residue_atom_id}{'label_alt_id'} = $alt_group_id;
        $residue_site->{$residue_atom_id}{'connections'} =
            $atom_site->{$residue_origin_atom_id}{'connections'};
        $atom_site->{$residue_origin_atom_id} =
            $residue_site->{$residue_atom_id};
    }

    return;
}

#
# Adds hydrogens to the molecule.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry. For example, hydrogens of
#     methyl group that has only one known connection are not added, but with two
#     connections - are;
#     $options->{use_existing_connections} - does not overwrite atom connections;
#     $options->{use_existing_hybridizations} - does not overwrite atom
#     hybridizations;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections;
#     $options->{exclude_by_atom_name} - list of atom names that are excluded
#     from being added hydrogens to;
#     $options->{alt_group_id} - alternative group id that is used to distinguish
#     pseudo atoms.
# Output:
#     %hydrogen_site - atom data structure with added hydrogens.
#

sub add_hydrogens
{
    my ( $atom_site, $options ) = @_;

    my ( $add_only_clear_positions, $use_existing_connections,
         $use_existing_hybridizations, $reference_atom_site,
         $exclude_by_atom_name, $exclude_by_atom_ids,
         $last_atom_id, $alt_group_id,
         $use_origins_alt_group_id ) =
        ( $options->{'add_only_clear_positions'},
          $options->{'use_existing_connections'},
          $options->{'use_existing_hybridizations'},
          $options->{'reference_atom_site'},
          $options->{'exclude_by_atom_name'},
          $options->{'exclude_by_atom_ids'},
          $options->{'last_atom_id'},
          $options->{'alt_group_id'},
          $options->{'use_origins_alt_group_id'}, );

    $add_only_clear_positions //= 0;
    $use_existing_connections //= 0;
    $use_existing_hybridizations //= 0;
    $reference_atom_site //= $atom_site;
    $exclude_by_atom_name //= [];
    $exclude_by_atom_ids //= [];
    $alt_group_id //= '1';
    $use_origins_alt_group_id //= 0;

    my %atom_site = %{ $atom_site };

    if( ! $use_existing_connections ) { connect_atoms( \%atom_site ) };

    if( ! $use_existing_hybridizations ) { hybridization( \%atom_site ) };

    my %hydrogen_site;
    $last_atom_id //= max( keys %{ $atom_site } );

    for my $atom_id ( sort { $a <=> $b } keys %atom_site ) {
        my $atom_name = $atom_site{$atom_id}{'label_atom_id'};

        next if any { $_ eq $atom_id } @{ $exclude_by_atom_ids };
        next if any { $_ eq $atom_name } @{ $exclude_by_atom_name };

        my $residue_name = $atom_site{$atom_id}{'label_comp_id'};

        my $hydrogen_names = $General::HYDROGEN_NAMES{$residue_name}{$atom_name};

        if( ! $hydrogen_names ) { next; }; # Exits early if there should be no
                                           # hydrogens connected to the atom.

        my $hybridization = $atom_site->{$atom_id}{'hybridization'};

        # Decides how many and what hydrogens should be added according to the
        # quantity of bonds and hydrogen atoms that should be connected to the
        # target atom.
        my @connection_ids = ();
        if( exists $reference_atom_site->{"$atom_id"}{'connections'} ) {
            @connection_ids =
                @{ $reference_atom_site->{"$atom_id"}{'connections'} };
        }

        # TODO: should be pre-determined as constant variable.
        my @mandatory_residue_atoms =
            @{ $General::RESIDUE_ATOMS{$residue_name}{'mandatory'} };
        my @mandatory_connections = ();
        for my $mandatory_atom ( @mandatory_residue_atoms ) {
            if( any { $mandatory_atom eq $_ }
                   @{ $General::CONNECTIVITY{$residue_name}{$atom_name} } ) {
                push @mandatory_connections, $mandatory_atom;
            }
        }

        # Hydrogens cannot be added if there is a missing information about
        # mandatory atom connections.
        next if( scalar @connection_ids < scalar @mandatory_connections );

        my @connection_names =
            map { $reference_atom_site->{"$_"}{'label_atom_id'} }
                @connection_ids;
        my @missing_hydrogens;
        for my $hydrogen_name ( @{ $hydrogen_names } ) {
            if( ! any { /$hydrogen_name/sxm } @connection_names ) {
                push @missing_hydrogens, $hydrogen_name;
            }
        }

        # Exits early if there are no spare hydrogens to add.
        if( scalar @missing_hydrogens == 0 ) { next; };

        #            sp3                       sp2               sp
        #
        #            Up(2)                     Up(2)             Up(2)
        # z          |                         |                 |
        # |_y      Middle(1) __ Right(3)     Middle(1)         Middle(1)
        # /         / \                       / \                |
        # x    Left(4) Back(5)           Left(4) Right(3)        Down(3)
        #
        # Depending on hybridization and present bond connections, adds missing
        # hydrogens.
        my %hydrogen_coord = map { $_ => undef } @missing_hydrogens;

        if( $hybridization eq 'sp3' ) {
            add_hydrogens_sp3( $atom_site, $atom_id, \%hydrogen_coord,
                               \@missing_hydrogens, $options );
        } elsif( $hybridization eq 'sp2' ) {
            add_hydrogens_sp2( $atom_site, $atom_id, \%hydrogen_coord,
                               \@missing_hydrogens, $options );
        } elsif( $hybridization eq 'sp' ) {
            add_hydrogens_sp( $atom_site, $atom_id, \%hydrogen_coord,
                              \@missing_hydrogens, $options );
        }

        # Each coordinate of atoms is transformed by transformation
        # matrix and added to %hydrogen_site.
        for my $hydrogen_name ( sort { $a cmp $b } keys %hydrogen_coord ) {
            if( $hydrogen_coord{$hydrogen_name} ) {
            # Adds necessary PDBx entries to pseudo atom site.
            $last_atom_id++;
            create_pdbx_entry(
                { 'atom_site' => \%hydrogen_site,
                  'id' => $last_atom_id,
                  'type_symbol' => 'H',
                  'label_atom_id' => $hydrogen_name,
                  'label_alt_id' =>
                      $use_origins_alt_group_id ?
                      $atom_site->{$atom_id}{'label_alt_id'} : $alt_group_id,
                  'label_comp_id' => $atom_site->{$atom_id}{'label_comp_id'},
                  'label_asym_id' => $atom_site->{$atom_id}{'label_asym_id'},
                  'label_entity_id' => $atom_site->{$atom_id}{'label_entity_id'},
                  'label_seq_id' => $atom_site{$atom_id}{'label_seq_id'},
                  'cartn_x' =>
                      sprintf( $SIG_FIGS_MIN,
                               $hydrogen_coord{$hydrogen_name}->[0][0] ),
                  'cartn_y' =>
                      sprintf( $SIG_FIGS_MIN,
                               $hydrogen_coord{$hydrogen_name}->[1][0] ),
                  'cartn_z' =>
                      sprintf( $SIG_FIGS_MIN,
                               $hydrogen_coord{$hydrogen_name}->[2][0] ),
                  'pdbx_PDB_model_num' =>
                      $atom_site->{$atom_id}{'pdbx_PDB_model_num'},
                } );
            # Adds additional pseudo-atom flag for future filtering.
            $hydrogen_site{$last_atom_id}{'is_pseudo_atom'} = 1;
            # Adds atom id that pseudo atoms was made of.
            $hydrogen_site{$last_atom_id}{'origin_atom_id'} = $atom_id;
            # Marks origin atom id as connection.
            $hydrogen_site{$last_atom_id}{'connections'} = [ $atom_id ];
            }
        }
    }

    return \%hydrogen_site;
}

#
# Adds hydrogens to the atoms which are sp3 hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp3
{
    my ( $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;

    my ( $add_only_clear_positions, $reference_atom_site ) = (
        $options->{'add_only_clear_positions'},
        $options->{'reference_atom_site'},
    );

    $add_only_clear_positions //= 0;
    $reference_atom_site //= $atom_site;

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $General::ATOMS{$atom_type}{'covalent_radius'}{'length'}[0] +
        $General::ATOMS{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids = @{ $reference_atom_site->{"$atom_id"}{'connections'} };
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my $lone_pair_count = $General::ATOMS{$atom_type}{'lone_pairs'};

    if( scalar( @connection_ids ) == 3 ) {
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord,
             $right_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]},
                   $atom_coord{$connection_ids[2]}, );

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
        #        alpha = arccos( ( - 4 - 2 * cos( beta ) -
        #                                2 * cos( gamma ) -
        #                                2 * cos( delta ) ) / 6 )
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
            bond_angle( [ $right_atom_coord,
                          $mid_atom_coord,
                          $left_atom_coord ] );

        # Calculates what common angle between hydrogen and rest of the
        # atoms should be.
        my $hydrogen_angle =
            acos( ( - 4 -
                    2 * cos( $up_mid_right_angle ) -
                    2 * cos( $up_mid_left_angle ) -
                    2 * cos( $right_mid_left_angle) ) /
                  6 );

        # Determines dihedral angle between left and right atoms. Then
        # splits rest of the 2 * pi angle into two equal parts.
        my $dihedral_angle =
            dihedral_angle( [ $left_atom_coord,
                              $up_atom_coord,
                              $mid_atom_coord,
                              $right_atom_coord ] );
        if( abs( $dihedral_angle ) < ( 3 * $PI / 4 ) ) {
            if( $dihedral_angle < 0 ) {
                $dihedral_angle = ( 2 * $PI + $dihedral_angle ) / 2;
            } else {
                $dihedral_angle = - ( 2 * $PI - $dihedral_angle ) / 2;
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
                   'global' ) };

        ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
            @{ mult_matrix_product(
                   [ $transf_matrix,
                     [ [ $bond_length *
                         cos( $PI / 2 - $dihedral_angle ) *
                         sin $hydrogen_angle ],
                       [ $bond_length *
                         sin( $PI / 2 - $dihedral_angle ) *
                         sin $hydrogen_angle ],
                       [ $bond_length *
                         cos $hydrogen_angle ],
                       [ 1 ] ] ] ) };
    } elsif( scalar( @connection_ids ) == 2 &&
           ! ( $add_only_clear_positions && $lone_pair_count > 0 ) ) {
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]}, );

        my $bond_angle =
            bond_angle( [ $up_atom_coord,
                          $mid_atom_coord,
                          $left_atom_coord ] );

        # This time is only one defined angle. And the calculation
        # changes:
        #
        #       acos( ( - 4 - 2 * cos( $bond_angle ) ) / 10 )
        #
        my $hydrogen_angle = acos( ( - 4 - 2 * cos $bond_angle ) / 10 );

        # Generates transformation matrix for transfering atoms to local
        # reference frame.
        my ( $transf_matrix ) =
            @{ switch_ref_frame( $mid_atom_coord,
                                 $up_atom_coord,
                                 $left_atom_coord,
                                 'global' ) };

        # Adds hydrogen first to both atoms that have 0 or 1 electron
        # pairs.
        if( scalar @{ $missing_hydrogens } >= 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 7 * $PI / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             sin( 7 * $PI / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             cos $hydrogen_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        # Additional hydrogen is added only to the atom that has no
        # electron pairs.
        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( - 1 * $PI / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             sin( - 1 * $PI / 6 ) *
                             sin $hydrogen_angle ],
                           [ $bond_length *
                             cos $hydrogen_angle ],
                           [ 1 ] ] ] ) };
        }

    } elsif( scalar @connection_ids == 1 && ! $add_only_clear_positions ){
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $side_coord ) = # Coordinate that only will be used for
                 # defining a local reference frame.
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   [ $atom_site->{$atom_id}{'Cartn_x'},
                     $atom_site->{$atom_id}{'Cartn_y'} + 1,
                     $atom_site->{$atom_id}{'Cartn_z'}, ], );

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
            $bond_angle = ( 109.5 - $lone_pair_count * 2.5 ) * $PI / 180;
        } else {
            $bond_angle = 109.5 * $PI / 180;
        }

        # Adds hydrogens according to the quantity of lone pairs.
        if( scalar @{ $missing_hydrogens } >= 3 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length * sin $bond_angle ],
                           [ 0 ],
                           [ $bond_length * cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } >= 2 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 2 * $PI / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 2 * $PI / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 4 * $PI / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             sin( 4 * $PI / 3 ) *
                             sin $bond_angle ],
                           [ $bond_length *
                             cos $bond_angle ],
                           [ 1 ] ] ] ) };
        }
    }

    return;
}

#
# Adds hydrogens to the atoms which are sp2 hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp2
{
    my ( $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;

    my ( $add_only_clear_positions, $reference_atom_site ) = (
        $options->{'add_only_clear_positions'},
        $options->{'reference_atom_site'}
    );

    $add_only_clear_positions //= 0;
    $reference_atom_site //= $atom_site;

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $General::ATOMS{$atom_type}{'covalent_radius'}{'length'}[1] +
        $General::ATOMS{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids = @{ $reference_atom_site->{"$atom_id"}{'connections'} };
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my $lone_pair_count = $General::ATOMS{$atom_type}{'lone_pairs'};

    # Depending on quantity of atoms connections, adds hydrogens.
    if( scalar @connection_ids == 2 ) {
        # Calculates current angle between atoms that are connected to
        # target atom.
        my ( $up_atom_coord,
             $mid_atom_coord,
             $left_atom_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   $atom_coord{$connection_ids[1]}, );

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
        #                  angle = ( 360 - beta ) / 2
        #
        # where beta is angle between two given bonds.
        my $bond_angle =
            ( 2 * $PI - bond_angle( [ $up_atom_coord,
                                      $mid_atom_coord,
                                      $left_atom_coord ] ) ) / 2;

        # Hydrogen is placed by placing hydrogen colinearly and then
        # rotating according to bond angle.
        ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
            @{ mult_matrix_product(
                   [ $transf_matrix,
                     [ [ $bond_length *
                         cos( - 0.5 * $PI ) *
                         sin $bond_angle ],
                       [ $bond_length *
                         sin( - 0.5 * $PI ) *
                         sin $bond_angle ],
                       [ $bond_length *
                         cos $bond_angle ],
                       [ 1 ] ] ] ) };

    } elsif( ( scalar @connection_ids == 1 ) &&
           ! ( $add_only_clear_positions && $lone_pair_count > 1 ) ) {
        my ( $up_atom_coord,
             $mid_atom_coord,
             $side_coord ) =
                 ( $atom_coord{$connection_ids[0]},
                   $atom_coord{$atom_id},
                   [ $atom_site->{$atom_id}{'Cartn_x'},
                     $atom_site->{$atom_id}{'Cartn_y'} + 1,
                     $atom_site->{$atom_id}{'Cartn_z'}, ], );

        # If terminal atom belongs to conjugated system, hydrogens are
        # added not to violate rule where atoms should be in one plain.
        my @second_neighbours =
            grep { ! /$atom_id/smx }
            map { @{ $reference_atom_site->{$_}{'connections'} } }
                @connection_ids;

        for my $second_neighbour ( @second_neighbours ) {
            my $second_hybridization =
                $reference_atom_site->{$second_neighbour}{'hybridization'};
            if( $second_hybridization eq 'sp2' ||
                $second_hybridization eq 'sp' ) {
                $side_coord =
                    [ $reference_atom_site->{$second_neighbour}{'Cartn_x'},
                      $reference_atom_site->{$second_neighbour}{'Cartn_y'},
                      $reference_atom_site->{$second_neighbour}{'Cartn_z'}, ];
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

        if( scalar @{ $missing_hydrogens } >= 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( -0.5 * $PI ) *
                             sin( 120 * $PI / 180 ) ],
                           [ $bond_length *
                             sin( -0.5 * $PI ) *
                             sin( 120 * $PI / 180 ) ],
                           [ $bond_length *
                             cos( 120 * $PI / 180 ) ],
                           [ 1 ] ] ] ) };
            shift @{ $missing_hydrogens };
        }

        if( scalar @{ $missing_hydrogens } == 1 ) {
            ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
                @{ mult_matrix_product(
                       [ $transf_matrix,
                         [ [ $bond_length *
                             cos( 0.5 * $PI ) *
                             sin( 120 * $PI / 180 ) ],
                           [ $bond_length *
                             sin( 0.5 * $PI ) *
                             sin( 120 * $PI / 180 ) ],
                           [ $bond_length *
                             cos( 120 * $PI / 180 ) ],
                           [ 1 ] ] ] ) };
        }
    }

    return;
}

#
# Adds hydrogens to the atoms which are sp hybridized.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $atom_id - atom id;
#     $hydrogen_coord - hydrogen coordinates from previous calculations - hash
#     of arrays;
#     $missing_hydrogens - list of hydrogens that are missing;
#     $options->{add_only_clear_postions} - adds only those atoms that have clear
#     positions that can be determined by geometry;
#     $options->{reference_atom_site} - atom site data structure that is used
#     for determining atom connections.
# Outputs:
#     appends hydrogen coordinates to $hydrogen_coord variable.
#

sub add_hydrogens_sp
{
    my ( $atom_site, $atom_id, $hydrogen_coord, $missing_hydrogens,
         $options ) = @_;
    my ( $reference_atom_site ) = ( $options->{'reference_atom_site'} );

    $reference_atom_site //= $atom_site;

    my $atom_type = $atom_site->{$atom_id}{'type_symbol'};

    my $bond_length =
        $General::ATOMS{$atom_type}{'covalent_radius'}{'length'}[1] +
        $General::ATOMS{'H'}{'covalent_radius'}{'length'}[0];

    my @connection_ids = @{ $reference_atom_site->{"$atom_id"}{'connections'} };
    my %atom_coord =
        %{ filter( { 'atom_site' => $reference_atom_site,
                     'include' => { 'id' => \@connection_ids },
                     'data' => [ 'Cartn_x', 'Cartn_y', 'Cartn_z' ],
                     'data_with_id' => 1 } ) };

    my ( $up_atom_coord,
         $mid_atom_coord,
         $side_coord ) =
             ( $atom_coord{$connection_ids[0]},
               $atom_coord{$atom_id},
               [ $atom_site->{$atom_id}{'Cartn_x'},
                 $atom_site->{$atom_id}{'Cartn_y'} + 1,
                 $atom_site->{$atom_id}{'Cartn_z'}, ], );

    # Generates transformation matrix for transfering atoms to local
    # reference frame.
    my ( $transf_matrix ) =
        @{ switch_ref_frame( $mid_atom_coord,
                             $up_atom_coord,
                             $side_coord,
                             'global' ) };

    ( $hydrogen_coord->{$missing_hydrogens->[0]} ) =
        @{ matrix_product(
               [ $transf_matrix,
                 [ [ 0 ],
                   [ 0 ],
                   [ - $bond_length ],
                   [ 1 ] ] ] ) };

    return;
}

1;
