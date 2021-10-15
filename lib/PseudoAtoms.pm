package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( calc_favourable_angle
                     calc_favourable_angles
                     calc_full_atom_energy
                     generate_library
                     generate_pseudo
                     generate_rotamer
                     library_to_csv
                     lowest_energy_state
                     replace_with_rotamer );

use B qw( svref_2object );
use Carp;
use Clone qw( clone );
use List::Util qw( max
                   shuffle );
use List::MoreUtils qw( any
                        uniq );
use Logging qw( info );
use threads;

use Combinatorics qw( permutation );
use ConnectAtoms qw( append_connections
                     connect_atoms
                     is_neighbour
                     is_second_neighbour
                     retains_connections );
use ForceField::Parameters;
use ForceField::Bonded qw( general );
use ForceField::NonBonded qw( general
                              hard_sphere
                              soft_sphere );
use Grid qw( grid_box
             identify_neighbour_cells );
use LinearAlgebra qw( mult_matrix_product );
use Measure qw( all_dihedral
                rmsd_sidechains );
use BondProperties qw( hybridization
                       rotatable_bonds
                       unique_rotatables );
use Multiprocessing qw( threading );
use PDBxParser qw( create_pdbx_entry
                   determine_residue_keys
                   filter_new
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
#     $args->{selection_state} - adds/changes '[local]_selection_state' to
#     specified value.
# Output:
#     $pseudo_atom_site - atom site data structure for pseudo-atoms with
#     additional 'conformation' attribute.
#

sub generate_pseudo
{
    my ( $args ) = @_;
    my ( $parameters, $atom_site, $atom_specifier, $angle_values, $last_atom_id,
         $alt_group_id, $selection_state ) =
        ( $args->{'parameters'}, $args->{'atom_site'}, $args->{'atom_specifier'},
          $args->{'angle_values'}, $args->{'last_atom_id'},
          $args->{'alt_group_id'}, $args->{'selection_state'}, );

    $last_atom_id //= max( keys %{ $atom_site } );
    $alt_group_id //= 1;

    my $sig_figs_max = $parameters->{'_[local]_constants'}{'sig_figs_max'};

    my %atom_site = %{ clone( $atom_site ) };
    my %pseudo_atom_site;

    my @atom_ids =
        @{ filter_new( \%atom_site,
                   { 'include' => $atom_specifier,
                     'return_data' => 'id' } ) };

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
        my @angle_names = grep { exists $angles{$residue_unique_key}{$_} }
                          sort { $a cmp $b }
                          keys %{ $angle_values };
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
                  'cartn_x' => sprintf( $sig_figs_max, $transf_atom_coord->[0][0] ),
                  'cartn_y' => sprintf( $sig_figs_max, $transf_atom_coord->[1][0] ),
                  'cartn_z' => sprintf( $sig_figs_max, $transf_atom_coord->[2][0] ),
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
            # Adds selection state if it is defined.
            if( defined $selection_state ) {
                $pseudo_atom_site{$last_atom_id}{'[local]_selection_state'} =
                    $selection_state;
            }
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
    my ( $parameters, $atom_site, $angle_values, $last_atom_id, $alt_group_id,
         $set_missing_angles_to_zero, $keep_origin_id, $keep_origin_alt_id ) =
        ( $args->{'parameters'}, $args->{'atom_site'}, $args->{'angle_values'},
          $args->{'last_atom_id'}, $args->{'last_atom_id'},
          $args->{'set_missing_angles_to_zero'}, $args->{'keep_origin_id'},
          $args->{'keep_origin_alt_id'}, );

    $last_atom_id //= max( keys %{ $atom_site } );
    $alt_group_id //= 1;
    $set_missing_angles_to_zero //= 0;
    $keep_origin_id //= 0;
    $keep_origin_alt_id //= 0; # Has higher priority than $alt_group_id.

    my %atom_site = %{ clone( $atom_site ) };
    my %rotamer_atom_site;

    for my $residue_unique_key ( keys %{ $angle_values } ) {
        my ( undef, undef, undef, $residue_alt_id ) =
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
                      'parameters' => $parameters,
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
#     $args->{angles} - angle data structure by which rotation is made:
#     {
#       'chi1' => {
#           'angle_start' => 0.0,
#           'angle_step' => 36.0,
#           'angle_end' => 360.0,
#     };
#     $args->{rmsd} - include RMSD calculations;
#     $args->{conf_model} - possible sidechain movements described by sidechain
#     model functions in SidechainModels.pm;
#     $args->{interactions} - interaction models described by functions in
#     AtomInteractions.pm;
#     $args->{parameters} - parameters that are passed to interaction function;
#     atom to have in the rotamer according to potential function;
#     $args->{threads} - number of threads.
# Output:
#     %library_atom_site - atom site data structure with additional data.
#

sub generate_library
{
    my ( $args ) = @_;
    my $parameters = $args->{'parameters'};
    my $atom_site = $args->{'atom_site'};
    my $residue_unique_keys = $args->{'residue_unique_keys'};
    my $include_interactions = $args->{'include_interactions'};
    my $angles = $args->{'angles'};
    my $rmsd = $args->{'rmsd'};
    my $conf_model = $args->{'conf_model'};
    my $interactions = $args->{'interactions'};
    my $threads = $args->{'threads'};
    my $program_called_by = $args->{'program_called_by'};
    my $options = $args->{'options'};

    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};
    my $interaction_atom_names = $parameters->{'_[local]_interaction_atom_names'};
    my $cutoff_atom = $parameters->{'_[local]_constants'}{'cutoff_atom'};

    $conf_model //= 'rotation_only';
    $threads //= 1;
    $include_interactions //= { 'label_atom_id' => $interaction_atom_names };
    $options //= {};

    # Selection of potential function.
    my %potential_functions =
        ( 'composite'   => { 'non_bonded' => \&ForceField::NonBonded::general,
                             'bonded'     => \&ForceField::Bonded::general},
          'hard_sphere' => { 'non_bonded' => \&hard_sphere },
          'soft_sphere' => { 'non_bonded' => \&soft_sphere }, );

    my %rotamer_library;

    my $atom_site_groups =
        split_by( { 'atom_site' => $atom_site,
                    'attributes' => [ 'pdbx_PDB_model_num', 'label_alt_id' ],
                    'append_dot' => 1  } );

    for my $atom_site_identifier ( sort keys %{ $atom_site_groups } ) {
        my ( $pdbx_model_num, $alt_id ) = split /,/, $atom_site_identifier;
        my $current_atom_site =
            filter_new( $atom_site,
                    { 'include' =>
                          {'id' => $atom_site_groups->{$atom_site_identifier}}});

        connect_atoms( $parameters, $current_atom_site );
        hybridization( $parameters, $current_atom_site );

        # Generates conformational models before checking for
        # clashes/interactions for given residues.
        if( $conf_model eq 'rotation_only' ) {
            rotation_only( $parameters, $current_atom_site );
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
                filter_new( $residue_site,
                        { 'include' => { 'label_atom_id' => [ 'CA' ] },
                          'return_data' => 'id' } )->[0];

            if( ! defined $atom_ca_id ) { next; }

            push @target_ca_ids, $atom_ca_id;
        }

        # Creates the grid box that has edge length of sum of all bonds of the
        # longest side-chain branch in arginine. Length: 3 * (C-C) + (C-N) + 2
        # * (C=N) + (N-H).
        # TODO: should consider using shorter distances, because bonds have limits
        # on maximum bending and having shorter edge length reduces calculation
        # time.
        my ( $grid_box, $target_cell_idxs ) =
            grid_box( $parameters, $current_atom_site, $edge_length_interaction,
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
                    filter_new( $current_atom_site,
                            { 'include' =>
                                  { 'pdbx_PDB_model_num' => [ $pdbx_model_num ],
                                    'label_alt_id' => [ $alt_id, q{.} ],
                                    'label_seq_id' => [ $residue_id ],
                                    'label_asym_id' => [ $residue_chain ] } } );
                my $residue_unique_key =
                    determine_residue_keys( $residue_site,
                                            { 'exclude_dot' => 1 } )->[0];

                my %interaction_site =
                    %{ filter_new( $current_atom_site,
                               { 'include' =>
                                     { 'id' => $neighbour_cells->{$cell},
                                       %{ $include_interactions } } } ) };

                # First, checks angles by step-by-step adding atoms to sidechains.
                # This is called growing side chain.
                my %options = %{ $options };
                my @allowed_angles =
                    @{ calc_favourable_angles(
                           { 'parameters' => $parameters,
                             'atom_site' => $current_atom_site,
                             'residue_unique_key' => $residue_unique_key,
                             'interaction_site' => \%interaction_site,
                             'angles' => $angles,
                             'non_bonded_potential' =>
                                 $potential_functions{$interactions}{'non_bonded'},
                             'bonded_potential' =>
                                 $potential_functions{$interactions}{'bonded'},
                             'threads' => $threads,
                             'options' => $options } ) };

                next if ! @allowed_angles;

                # Then, re-checks if each atom of the rotamer obey energy
                # cutoffs.
                my ( $allowed_angles, $energy_sums, $rmsds ) =
                    @{ threading(
                           \&calc_full_atom_energy,
                           { 'parameters' => $parameters,
                             'atom_site' => $current_atom_site,
                             'residue_unique_key' => $residue_unique_key,
                             'interaction_site' => \%interaction_site,
                             'non_bonded_potential' =>
                                 $potential_functions{$interactions}{'non_bonded'},
                             'bonded_potential' =>
                                 $potential_functions{$interactions}{'bonded'},
                             ( $rmsd ? ( 'rmsd' => 1 ): ()  ),
                             'options' => $options },
                           [ \@allowed_angles ],
                           $threads ) };

                # my ( $allowed_angles, $energy_sums ) =
                #     @{ calc_full_atom_energy(
                #            { 'parameters' => $parameters,
                #              'atom_site' => $current_atom_site,
                #              'residue_unique_key' => $residue_unique_key,
                #              'interaction_site' => \%interaction_site,
                #              'non_bonded_potential' =>
                #                  $potential_functions{$interactions}{'non_bonded'},
                #              'bonded_potential' =>
                #                  $potential_functions{$interactions}{'bonded'},
                #              ( $rmsd ? ( 'rmsd' => 1 ): ()  ),
                #              'options' => $options },
                #            [ @allowed_angles ] ) };

                for( my $i = 0; $i <= $#{ $allowed_angles }; $i++  ) {
                    my %angles =
                        map { my $angle_id = $_ + 1;
                              ( "chi$angle_id" => $allowed_angles->[$i][$_] ) }
                            ( 0..$#{ $allowed_angles->[$i] } );
                    my $rotamer_energy_sum = $energy_sums->[$i];
                    if( defined $rotamer_energy_sum ) {
                        push @{ $rotamer_library{"$residue_unique_key"} },
                            { 'angles' => \%angles,
                              'potential' => $interactions,
                              'potential_energy_value' => $energy_sums->[$i],
                              ( $rmsd ? ( 'rmsd' => $rmsds->[$i][-1] ) : () ) };
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
#     $args->{angles} - angle data structure by which rotation is made.
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of non-bonded atoms;
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of bonded atoms;
#     $args->{parameters} - parameters that are passed to interaction function;
#     $args->{threads} - number of threads.
# Output:
#     @allowed_angles - list of groups of allowed angles.
#     Ex.: ( [ 0.00, 3.14 ], [ 3.14, 6.28 ] ).
#

sub calc_favourable_angles
{
    my ( $args ) = @_;

    my ( $parameters, $atom_site, $residue_unique_key, $interaction_site,
         $angles, $small_angle, $non_bonded_potential, $bonded_potential,
         $threads, $rand_count, $rand_seed, $program_called_by, $verbose ) = (
        $args->{'parameters'},
        $args->{'atom_site'},
        $args->{'residue_unique_key'},
        $args->{'interaction_site'},
        $args->{'angles'},
        $args->{'small_angle'},
        $args->{'non_bonded_potential'},
        $args->{'bonded_potential'},
        $args->{'threads'},
        $args->{'options'}{'rand_count'},
        $args->{'options'}{'rand_seed'},
        $args->{'options'}{'program_called_by'},
        $args->{'options'}{'verbose'},
    );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    # TODO: look how separate $angles and $small_angle influence on the function.
    $small_angle //= 0.1 * 2 * $pi;
    $rand_seed //= 23;

    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $residue_unique_key, 1 );

    my ( $any_key ) = keys %{ $residue_site };
    my $residue_name = $residue_site->{$any_key}{'label_comp_id'};

    my $rotatable_bonds = rotatable_bonds( $residue_site );
    if( ! %{ $rotatable_bonds } ) { return []; }

    # Goes through each atom in side chain and calculates interaction
    # potential with surrounding atoms. CA and CB are non-movable atoms
    # so, they are marked as starting atoms.
    my $ca_atom_id =
        filter_new( $residue_site,
                { 'include' => { 'label_atom_id' => [ 'CA' ] },
                  'return_data' => 'id' } )->[0];
    my $cb_atom_id =
        filter_new( $residue_site,
                { 'include' => { 'label_atom_id' => [ 'CB' ] },
                  'return_data' => 'id' } )->[0];

    my @visited_atom_ids = ( $ca_atom_id, $cb_atom_id );
    my @next_atom_ids =
        grep { $_ ne $ca_atom_id }
            @{ $residue_site->{$cb_atom_id}{'connections'} };

    my @allowed_angles;
    my @allowed_energies;
    while( scalar( @next_atom_ids ) != 0 ) {
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my @default_allowed_angles;
            # TODO: last angle should be sorted with <=> by first removing chi
            # prefix. It will be important if large quantity of dihedral angles
            # are analyzed.
            my ( $last_angle_name ) =
                sort { $b cmp $a } keys %{ $rotatable_bonds->{$atom_id} };

            if( exists $angles->{$residue_name}{$last_angle_name} ) {
                @default_allowed_angles =
                    map { [ $_ ] } @{ $angles->{$residue_name}{$last_angle_name} };
            } elsif( exists $angles->{$residue_name}{'*'} ) {
                @default_allowed_angles =
                    map { [ $_ ] } @{ $angles->{$residue_name}{'*'} };
            } elsif( exists $angles->{'*'}{$last_angle_name} ) {
                @default_allowed_angles =
                    map { [ $_ ] } @{ $angles->{'*'}{$last_angle_name} };
            } elsif( exists $angles->{'*'}{'*'} ) {
                if( defined $rand_count && defined $rand_seed ) {
                    if( $rand_count > scalar @{$angles->{'*'}{'*'}} ) {
                        die 'number of randomly selected angles is greater that ' .
                            "possible angles.\n";
                    }
                    my @shuffled_idxs = shuffle( 0..$#{$angles->{'*'}{'*'}} );
                    @default_allowed_angles =
                        map { [ $angles->{'*'}{'*'}[$_] ] }
                            @shuffled_idxs[0..$rand_count-1];
                } else {
                    @default_allowed_angles =
                        map { [ $_ ] } @{ $angles->{'*'}{'*'} };
                }
            } else {
                @default_allowed_angles =
                    map { [ $_ ] }
                       @{ sample_angles( $parameters,
                                         [ [ 0, 2 * $pi ] ],
                                         $small_angle ) };
            }

            my @default_allowed_energies = map { [ 0 ] } @default_allowed_angles;

            # Adds more angle combinations if there are more than one
            # rotatable bonds.
            if( @allowed_angles &&
                scalar( @{ $allowed_angles[0] } ) <
                scalar( keys %{ $rotatable_bonds->{$atom_id} } ) ) {
                @allowed_angles =
                    @{ permutation( 2, [], [ \@allowed_angles,
                                             \@default_allowed_angles ], [] ) };
                @allowed_energies =
                    @{ permutation( 2, [], [ \@allowed_energies,
                                             \@default_allowed_energies ], [] ) };
                # Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
                @allowed_angles =
                    map { [ @{ $_->[0] }, @{ $_->[1] } ] } @allowed_angles;
                @allowed_energies =
                    map { [ $_->[0][0] ] } @allowed_energies;
            } elsif( ! @allowed_angles ) {
                @allowed_angles = @default_allowed_angles;
                @allowed_energies = @default_allowed_energies;
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            # Marks neighbouring atoms.
            push @neighbour_atom_ids, @{ $atom_site->{$atom_id}{'connections'} };

            # Starts calculating potential energy.
            my ( $next_allowed_angles, $next_allowed_energies ) =
                @{ threading(
                       \&calc_favourable_angle,
                       { 'parameters' => $parameters,
                         'atom_site' => $atom_site,
                         'atom_id' => $atom_id,
                         'interaction_site' => $interaction_site,
                         'non_bonded_potential' => $non_bonded_potential,
                         'bonded_potential' => $bonded_potential },
                       [ \@allowed_angles, \@allowed_energies ],
                       $threads ) };

            # my ( $next_allowed_angles, $next_allowed_energies ) =
            #     @{ calc_favourable_angle(
            #            { 'parameters' => $parameters,
            #              'atom_site' => $atom_site,
            #              'atom_id' => $atom_id,
            #              'interaction_site' => $interaction_site,
            #              'non_bonded_potential' => $non_bonded_potential,
            #              'bonded_potential' => $bonded_potential },
            #            [ \@allowed_angles, \@allowed_energies ] ) };

            if( scalar @{ $next_allowed_angles } > 0 ) {
                @allowed_angles = @{ $next_allowed_angles };
                @allowed_energies = @{ $next_allowed_energies };
                print info(
                    { message =>
                          $residue_site->{$atom_id}{'pdbx_PDB_model_num'} . " " .
                          $residue_site->{$atom_id}{'label_asym_id'} . " " .
                          $residue_site->{$atom_id}{'label_seq_id'} . " " .
                          $residue_site->{$atom_id}{'label_alt_id'} . " " .
                          "${residue_name} " .
                          $residue_site->{$atom_id}{'label_atom_id'} . " " .
                          "${last_angle_name} " . scalar( @allowed_angles ) . "\n",
                      program => $program_called_by }
                    ) if $verbose;
            } else {
                return [];
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
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of non-bonded atoms;
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of bonded atoms;
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

    my ( $parameters, $atom_site, $atom_id, $interaction_site,
         $non_bonded_potential, $bonded_potential, $options ) = (
        $args->{'parameters'},
        $args->{'atom_site'},
        $args->{'atom_id'},
        $args->{'interaction_site'},
        $args->{'non_bonded_potential'},
        $args->{'bonded_potential'},
        $args->{'options'},
    );

    my $energy_cutoff_atom= $parameters->{'_[local]_force_field'}{'cutoff_atom'};

    my %options = defined $options ? %{ $options } : ();
    $options{'atom_site'} = $atom_site;

    my @allowed_angles;
    my @allowed_energies;
    for( my $i = 0; $i <= $#{ $array_blocks->[0] }; $i++ ) {
        my $angles = $array_blocks->[0][$i];
        my $energies = $array_blocks->[1][$i][0];
        my %angles =
            map { my $angle_id = $_ + 1; ( "chi$angle_id" => [ $angles->[$_] ] )}
                ( 0..$#{ $angles } );

        my $pseudo_atom_site =
            generate_pseudo( { 'parameters' => $parameters,
                               'atom_site' => $atom_site,
                               'atom_specifier' => { 'id' => [ "$atom_id" ] },
                               'angle_values' => \%angles } );
        my $pseudo_atom_id = ( keys %{ $pseudo_atom_site } )[0];
        my $pseudo_origin_id =
            $pseudo_atom_site->{$pseudo_atom_id}{'origin_atom_id'};

        my $potential_energy = 0; # TODO: look if here should be zeros.
        my $potential_sum = 0;

        # Calculation of potential energy of non-bonded atoms.
        foreach my $interaction_id ( keys %{ $interaction_site } ) {
            if( ( ! is_neighbour( $atom_site,
                                  $pseudo_origin_id,
                                  $interaction_id ) ) &&
                ( ! is_second_neighbour( $atom_site,
                                         $pseudo_origin_id,
                                         $interaction_id ) ) ) {
                $potential_energy = $non_bonded_potential->(
                    $parameters,
                    $pseudo_atom_site->{$pseudo_atom_id},
                    $atom_site->{$interaction_id},
                    \%options
                );
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
#     $args->{angles} - angle data structure by which rotation is made.
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of non-bonded atoms;
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

    my ( $parameters, $atom_site, $residue_unique_key, $interaction_site,
         $non_bonded_potential, $bonded_potential, $rmsd, $options ) = (
        $args->{'parameters'},
        $args->{'atom_site'},
        $args->{'residue_unique_key'},
        $args->{'interaction_site'},
        $args->{'non_bonded_potential'},
        $args->{'bonded_potential'},
        $args->{'rmsd'},
        $args->{'options'},
    );

    my $energy_cutoff_atom = $parameters->{'_[local]_force_field'}{'cutoff_atom'};
    my $interaction_atom_names = $parameters->{'_[local]_interaction_atom_names'};

    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $residue_unique_key, 1 );

    # Checks for inter-atom interactions and determines if energies
    # comply with cutoffs.
    my @checkable_angles = @{ $array_blocks->[0] };
    my @allowed_angles;
    my @energy_sums;
    my @rmsd_averages;

  ALLOWED_ANGLES:
    for( my $i = 0; $i <= $#checkable_angles; $i++ ) {
        my %angles =
            map { my $angle_id = $_ + 1;
                  ( "chi$angle_id" => $checkable_angles[$i][$_] ) }
                ( 0..$#{ $checkable_angles[$i] } );

        my %rotamer_site = %{ $residue_site };
        replace_with_rotamer( $parameters, \%rotamer_site, $residue_unique_key,
                              \%angles );

        # connect_atoms($parameters, \%rotamer_site, {'only_covalent_radii' => 1});

        # next if !retains_connections( $parameters,$residue_site,\%rotamer_site );

        my @rotamer_atom_ids =
            sort keys %{ filter_new( \%rotamer_site,
                                 { 'exclude' =>
                                   { 'label_atom_id' =>
                                         $interaction_atom_names } } ) };

        # HACK: make sure that $interaction_site atom ids are updated by
        # %rotamer_site.
        my %rotamer_interaction_site = ( %{ $interaction_site }, %rotamer_site );

        # HACK: should connect_atoms() be used here?
        # connect_atoms( \%rotamer_interaction_site );

        $options->{'atom_site'} = \%rotamer_interaction_site;

        my $rotamer_energy_sum = 0;

        for my $rotamer_atom_id ( @rotamer_atom_ids ) {
            # Calculation of potential energy of bonded atoms.
            if( defined $bonded_potential ) {
                $rotamer_energy_sum += $bonded_potential->(
                    $parameters,
                    $rotamer_interaction_site{$rotamer_atom_id},
                    $options
                );
            }

            # Calculation of potential energy of non-bonded atoms.
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
                        $non_bonded_potential->(
                            $parameters,
                            $rotamer_interaction_site{$rotamer_atom_id},
                            $rotamer_interaction_site{$neighbour_atom_id},
                            $options );

                    next ALLOWED_ANGLES
                        if $rotamer_atom_energy > $energy_cutoff_atom;

                    $rotamer_energy_sum += $rotamer_atom_energy;
                }
            }
        }

        push @allowed_angles, $checkable_angles[$i];
        push @energy_sums, $rotamer_energy_sum;

        if( defined $rmsd ) {
            push @rmsd_averages,
                map { [ $_->{'value'} ] }
                   @{ rmsd_sidechains( $parameters, $residue_site,\%rotamer_site,
                                       { 'average' => 1 } ) };
        }
    }

    return [ \@allowed_angles, \@energy_sums, \@rmsd_averages ] ;
}

#
# Calculates lowest possible energy state that atom would have if it interacted
# with all surrounding atoms with best possible distances. Although, it is
# very unrealistic for atom to be in such state, but very useful when using in
# determining energy cutoff value for dead-end elimination algorithm.
# Input:
#     $atom_i - atom;
#     $surrounding_atoms - list of surrounding atoms;
#     $non_bonded_potential - reference to the potential function that is used for
#     calculating energy of non-bonded atoms;
#     $bonded_potential - reference to the potential function that is used for
#     calculating energy of bonded atoms;
#     $parameters - potential function parameters.
# Output:
#     $lowest_energy_sum - energy value.
#

# FIXME: add bonded potential.
sub lowest_energy_state
{
    my ( $parameters, $atom_i, $surrounding_atoms, $non_bonded_potential,
         $options ) = @_;

    $options //= {};

    my $lowest_energy_sum = 0;
    for my $atom_j ( @{ $surrounding_atoms } ) {
        $lowest_energy_sum +=
            $non_bonded_potential->( $parameters, $atom_i, $atom_j,
                                   { %{ $options }, ( 'is_optimal' => 1 ) } );
    }

    return $lowest_energy_sum;
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
    my ( $parameters, $atom_site, $residue_unique_key, $angle_values ) = @_;

    my ( undef, undef, undef, $alt_group_id ) = split /,/, $residue_unique_key;
    my $residue_site =
        generate_rotamer( { 'parameters' => $parameters,
                            'atom_site' => $atom_site,
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

1;
