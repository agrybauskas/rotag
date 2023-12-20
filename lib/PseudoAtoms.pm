package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( assign_hetatoms_with_struct_conn
                     assign_hetatoms_no_struct_conn
                     calc_favourable_angle
                     calc_favourable_angles
                     calc_full_atom_energy
                     generate_library
                     generate_pseudo
                     generate_rotamer
                     lowest_energy_state
                     replace_with_rotamer );
}

use Carp;
use Clone qw( clone );
use List::Util qw( max
                   shuffle );
use List::MoreUtils qw( any
                        uniq );
use Logging qw( info
                warning );
use threads;

use BondParameters qw( collect_bond_lengths
                       collect_bond_angles
                       collect_dihedral_angles
                       bendable_angles
                       rotatable_bonds
                       stretchable_bonds );
use BondProperties qw( hybridization );
use Combinatorics qw( permutation );
use ConnectAtoms qw( connect_atoms
                     connect_atoms_explicitly
                     disconnect_atoms_explicitly
                     is_neighbour
                     is_second_neighbour );
use ForceField::Parameters;
use ForceField::Bonded qw( general );
use ForceField::NonBonded qw( general
                              hard_sphere
                              soft_sphere );
use Grid qw( grid_box
             identify_neighbour_cells );
use LinearAlgebra qw( mult_matrix_product );
use Measure qw( around_distance
                rmsd_sidechains );
use Moieties qw( missing_atom_names );
use Multiprocessing qw( threading );
use PDBxParser qw( create_pdbx_entry
                   determine_residue_keys
                   filter_new
                   filter_by_unique_residue_key
                   split_by
                   unique_residue_key );
use Sampling qw( sample_angles );
use SidechainModels qw( rotation_translation );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates hetero pseudo-atoms from atom models that are written in list of
# analytical matrices.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser). Must have
#     any sidechain model function applied to it (see SidechainModels.pm);
#     $args->{atom_specifier} - hash of hashes for selecting atoms by attributes
#     (see PDBxParser.pm);
#     $args->{angle_and_bond_values} - hash of arrays that describe possible
#     values of dihedral angles, bond lengths and angles:
#     Ex.: { 'chi1' => [ 0, 0.4, 1.5, 2.0 ],
#            'chi2' => [ 0, 2 ] };
#     $args->{last_atom_id} - last atom id for assigning new ids for pseud
#     atoms;
#     $args->{alt_group_id} - alternative group id that is used to distinguish
#     pseudo atoms. Very useful when generating rotamers.
#     $args->{selection_state} - adds/changes '[local]_selection_state' to
#     specified value.
# Output:
#     $pseudo_hetatom_site - atom site data structure for pseudo-atoms with
#     additional 'conformation' attribute.
#

sub generate_pseudo
{
    my ( $args ) = @_;
    my ( $parameters, $atom_site, $atom_specifier, $bond_parameter_values,
         $last_atom_id, $alt_group_id, $selection_state ) =
        ( $args->{'parameters'}, $args->{'atom_site'}, $args->{'atom_specifier'},
          $args->{'bond_parameter_values'}, $args->{'last_atom_id'},
          $args->{'alt_group_id'}, $args->{'selection_state'} );

    $last_atom_id //= max( keys %{ $atom_site } );
    $alt_group_id //= 1;

    my $sig_figs_max = $parameters->{'_[local]_constants'}{'sig_figs_max'};

    my %atom_site = %{ clone( $atom_site ) };
    my %pseudo_atom_site;

    my @atom_ids = @{ filter_new( \%atom_site,
                                  { 'include' => $atom_specifier,
                                    'return_data' => 'id' } ) };

    for my $atom_id ( @atom_ids ) {
        my $conformation = $atom_site{"$atom_id"}{'conformation'};

        confess "atom with id $atom_id lacks 'conformation' key."
            if ! defined $conformation;

        # Calculates current dihedral angles of rotatable bonds, bond lengths
        # and angle between bonds. Will be used for reseting dihedral angles to 0
        # degree angle and pinpointing the correct bond length and angle changes.
        my $residue_unique_key = unique_residue_key( $atom_site{$atom_id} );

        my %bond_parameters =
            ( %{ collect_dihedral_angles( { $atom_id => $atom_site->{$atom_id} } )
                 ->{$residue_unique_key} },
              %{ collect_bond_lengths( { $atom_id => $atom_site->{$atom_id} } )
                 ->{$residue_unique_key} },
              %{ collect_bond_angles( { $atom_id => $atom_site->{$atom_id} } )
                 ->{$residue_unique_key} } );

        # Adjust changes to the existing values of the bond and angle parameters.
        my @bond_parameter_names = sort keys %bond_parameters;
        my @bond_parameter_values = ();
        for my $bond_parameter_name ( @bond_parameter_names  ) {
            if( exists $bond_parameter_values->{"$bond_parameter_name"} ) {
                push @bond_parameter_values,
                    [ map { $_ - $bond_parameters{"$bond_parameter_name"}{'value'} }
                          @{ $bond_parameter_values->{"$bond_parameter_name"} } ];
            } else {
                push @bond_parameter_values, [ 0.0 ];
            }
        }

        # Iterates through combinations of angles, lengths and evaluates
        # conformational model.
        for my $bond_parameter_comb (
            @{ permutation( scalar( @bond_parameter_names ), [],
                            \@bond_parameter_values, [] ) } ){

            my %bond_parameter_values =
                map {( $bond_parameter_names[$_] => $bond_parameter_comb->[$_] )}
                    0..$#bond_parameter_names;

            # Evaluates matrices.
            my ( $transf_atom_coord ) = @{
                mult_matrix_product( $conformation,
                                     { %bond_parameter_values, 'eta' => 0.0 } )
            };

            # Adds necessary PDBx entries to pseudo atom site.
            $last_atom_id++;
            create_pdbx_entry(
                { 'group_PDB' => $atom_site->{$atom_id}{'group_PDB'},
                  'atom_site' => \%pseudo_atom_site,
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
            $pseudo_atom_site{$last_atom_id}{'connections_hetatom'} =
                $atom_site{$atom_id}{'connections_hetatom'};
            $pseudo_atom_site{$last_atom_id}{'conformation'} =
                $atom_site{$atom_id}{'conformation'};

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

        for my $atom_id ( sort { $a <=> $b } keys %{ $residue_site } ) {
            my $rotatable_bonds = $atom_site->{$atom_id}{'rotatable_bonds'};
            my $stretchable_bonds = $atom_site->{$atom_id}{'stretchable_bonds'};
            my $bendable_angles = $atom_site->{$atom_id}{'bendable_angles'};

            my %bond_parameters = (
                ( defined $rotatable_bonds ? %{ $rotatable_bonds } : () ),
                ( defined $stretchable_bonds ? %{ $stretchable_bonds } : () ),
                ( defined $bendable_angles ? %{ $bendable_angles } : () )
            );

            if( ! %bond_parameters ) { next; }

            my %angles;
            for my $angle_name ( keys %bond_parameters ) {
                if( exists $angle_values->{"$residue_unique_key"}{$angle_name} &&
                    defined $angle_values->{"$residue_unique_key"}{$angle_name}){
                    $angles{$angle_name} =
                        [ $angle_values->{"$residue_unique_key"}{$angle_name} ];
                } else {
                    if( $set_missing_angles_to_zero ) {
                        $angles{$angle_name} = [ 0.0 ];
                    } else {
                        $angles{$angle_name} =
                            [ $bond_parameters{$angle_name}{'value'} ];
                    }
                }
            }

            %rotamer_atom_site =
                ( %rotamer_atom_site,
                  %{ generate_pseudo( {
                      'parameters' => $parameters,
                      'atom_site' => { ( %atom_site, %rotamer_atom_site ) },
                      'atom_specifier' => { 'id' => [ $atom_id ] },
                      'bond_parameter_values' => \%angles,
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
#     $args->{bond_parameters} - bond parameter data structure by which rotation
#     or/and translation is made:
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
    my $struct_conn = $args->{'struct_conn'};
    my $residue_unique_keys = $args->{'residue_unique_keys'};
    my $include_interactions = $args->{'include_interactions'};
    my $include_hetatoms = $args->{'include_hetatoms'};
    my $bond_parameters = $args->{'bond_parameters'};
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

    $struct_conn //= {};
    $conf_model //= 'rotation_only';
    $threads //= 1;
    $include_interactions //= { 'label_atom_id' => $interaction_atom_names };
    $include_hetatoms //= 0;
    $options //= {};

    my $do_bond_torsion =
        is_bond_parameter_present( $bond_parameters, 'dihedral_angle' );
    my $do_bond_stretching = 0;
        is_bond_parameter_present( $bond_parameters, 'bond_length' );
    my $do_angle_bending = 0;
        is_bond_parameter_present( $bond_parameters, 'bond_angle' );

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

        if( $include_hetatoms ) {
            if( %{ $struct_conn } ) {
                assign_hetatoms_with_struct_conn(
                    $parameters, $current_atom_site, $struct_conn
                );
            } else {
                assign_hetatoms_no_struct_conn(
                    $parameters, $current_atom_site
                );
            }
        }

        if( $do_bond_torsion ) {
            rotatable_bonds( $parameters,
                             $current_atom_site,
                             { 'include_hetatoms' => $include_hetatoms } );
        }
        if( $do_bond_stretching ) {
            stretchable_bonds( $parameters,
                               $current_atom_site,
                               { 'include_hetatoms' => $include_hetatoms } );
        }
        if( $do_angle_bending ) {
            bendable_angles( $parameters,
                             $current_atom_site,
                             { 'include_hetatoms' => $include_hetatoms } );
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
                my $residue_name =
                    $current_atom_site->{$ca_atom_id}{'label_comp_id'};
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

                my $dihedral_angles =
                    collect_dihedral_angles( $residue_site )->{$residue_unique_key};
                my $bendable_angles =
                    collect_bond_angles( $residue_site )->{$residue_unique_key};
                my $stretchable_bonds =
                    collect_bond_lengths( $residue_site )->{$residue_unique_key};

                my %bond_parameters = (
                    ( defined $dihedral_angles ? %{ $dihedral_angles } : () ),
                    ( defined $stretchable_bonds ? %{ $stretchable_bonds } : () ),
                    ( defined $bendable_angles ? %{ $bendable_angles } : () ),
                );

                my @bond_parameter_names =
                    sort { $bond_parameters{$a}{'order'} <=>
                           $bond_parameters{$b}{'order'} }
                    keys %bond_parameters;

                my @missing_atom_names =
                    @{ missing_atom_names( $parameters, $residue_site ) };

                if( @missing_atom_names ) {
                    warn 'missing ' .
                        ( $#missing_atom_names > 0 ? 'atoms' : 'atom' ) .
                        ' in the residue ' . uc( $residue_name ) . $residue_id .
                        ': ' . join( ', ',  @missing_atom_names ) . "\n";
                }

                # Generates conformational models before checking for
                # clashes/interactions for given residues.
                if( $conf_model eq 'rotation_only' ) {
                    rotation_translation( $parameters, $residue_site );
                } else {
                    confess 'conformational model was not defined.';
                }

                my %interaction_site =
                    %{ filter_new( $current_atom_site,
                               { 'include' =>
                                     { 'id' => $neighbour_cells->{$cell},
                                       %{ $include_interactions } } } ) };

                # First, checks bond, dihedral angles and bond length changes
                # by step-by-step adding atoms to sidechains. This is called
                # growing side chain.
                my %options = %{ $options };
                my @allowed_angles =
                    @{ calc_favourable_angles(
                           { 'parameters' => $parameters,
                             'atom_site' => $current_atom_site,
                             'residue_unique_key' => $residue_unique_key,
                             'interaction_site' => \%interaction_site,
                             'bond_parameters' => $bond_parameters,
                             'include_hetatoms' => $include_hetatoms,
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
                              ( $bond_parameter_names[$_] => $allowed_angles->[$i][$_] ) }
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
# Calculates favourable rotamer dihedral, bond angles and bond length changes
# for a given residue.
# Input:
#     $args->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $args->{residue_unique_key} - unique residue key
#     (see PDBxParser::unique_residue_key);
#     $args->{interaction_site} - atom data structure that is included into
#     energy calculations;
#     $args->{bond_parameters_angles} - bond parameter data structure by which
#     conformation is made.
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of non-bonded atoms;
#     $args->{non_bonded_potential} - reference to the potential function that is
#     used for calculating energy of bonded atoms;
#     $args->{parameters} - parameters that are passed to interaction function;
#     $args->{threads} - number of threads.
# Output:
#     @allowed_bond_parameters - list of groups of allowed bond parameters.
#     Ex.: ( [ 0.00, 3.14 ], [ 3.14, 6.28 ] ).
#

sub calc_favourable_angles
{
    my ( $args ) = @_;

    my ( $parameters, $atom_site, $residue_unique_key, $interaction_site,
         $bond_parameters, $include_hetatoms, $bond_parameter_count,
         $non_bonded_potential, $bonded_potential, $threads, $rand_count,
         $rand_seed, $program_called_by, $verbose ) = (
        $args->{'parameters'},
        $args->{'atom_site'},
        $args->{'residue_unique_key'},
        $args->{'interaction_site'},
        $args->{'bond_parameters'},
        $args->{'include_hetatoms'},
        $args->{'bond_parameter_count'},
        $args->{'non_bonded_potential'},
        $args->{'bonded_potential'},
        $args->{'threads'},
        $args->{'options'}{'rand_count'},
        $args->{'options'}{'rand_seed'},
        $args->{'options'}{'program_called_by'},
        $args->{'options'}{'verbose'},
    );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    $bond_parameter_count //= 20;
    $rand_seed //= 23;
    $include_hetatoms //= 0;

    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $residue_unique_key, 1 );

    my ( $any_key ) = keys %{ $residue_site };
    my $residue_name = $residue_site->{$any_key}{'label_comp_id'};

    my $dihedral_angles =
        collect_dihedral_angles( $residue_site )->{$residue_unique_key};
    my $bendable_angles =
        collect_bond_angles( $residue_site )->{$residue_unique_key};
    my $stretchable_bonds =
        collect_bond_lengths( $residue_site )->{$residue_unique_key};

    my %bond_parameters = (
        ( defined $dihedral_angles ? %{ $dihedral_angles } : () ),
        ( defined $stretchable_bonds ? %{ $stretchable_bonds } : () ),
        ( defined $bendable_angles ? %{ $bendable_angles } : () ),
    );

    if( ! %bond_parameters ) { return []; }

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
             ( @{ $residue_site->{$cb_atom_id}{'connections'} },
               @{ $include_hetatoms &&
                  defined $residue_site->{$cb_atom_id}{'connections_hetatom'} ?
                  $residue_site->{$cb_atom_id}{'connections_hetatom'} :
                  [] } );

    my @allowed_bond_parameters;
    my @allowed_energies;
    while( scalar( @next_atom_ids ) != 0 ) {
        my @neighbour_atom_ids;
        for my $atom_id ( @next_atom_ids ) {
            my @default_allowed_bond_parameters;

            my ( $last_angle_name ) =
                sort { $bond_parameters{$b}{'order'} <=>
                       $bond_parameters{$a}{'order'} }
                keys %bond_parameters;

            if( exists $bond_parameters->{$residue_name}{$last_angle_name} ) {
                @default_allowed_bond_parameters =
                    map { [ $_ ] }
                       @{ $bond_parameters->{$residue_name}{$last_angle_name} };
            } elsif( exists $bond_parameters->{$residue_name}{'*'} ) {
                @default_allowed_bond_parameters =
                    map { [ $_ ] }
                       @{ $bond_parameters->{$residue_name}{'*'} };
            } elsif( exists $bond_parameters->{'*'}{$last_angle_name} ) {
                @default_allowed_bond_parameters =
                    map { [ $_ ] }
                       @{ $bond_parameters->{'*'}{$last_angle_name} };
            } elsif( exists $bond_parameters->{'*'}{'*'} ) {
                if( defined $rand_count && defined $rand_seed ) {
                    if( $rand_count > scalar @{$bond_parameters->{'*'}{'*'}} ) {
                        die 'number of randomly selected angles is greater ' .
                            "that possible angles.\n";
                    }
                    my @shuffled_idxs =
                        shuffle( 0..$#{$bond_parameters->{'*'}{'*'}} );
                    @default_allowed_bond_parameters =
                        map { [ $bond_parameters->{'*'}{'*'}[$_] ] }
                            @shuffled_idxs[0..$rand_count-1];
                } else {
                    @default_allowed_bond_parameters =
                        map { [ $_ ] } @{ $bond_parameters->{'*'}{'*'} };
                }
            } else {
                @default_allowed_bond_parameters =
                    map { [ $_ ] }
                       @{ sample_angles( [ [ 0, 2 * $pi ] ],
                                         $bond_parameter_count ) };
            }

            my @default_allowed_energies =
                map { [ 0 ] } @default_allowed_bond_parameters;

            # Adds more angle combinations if there are more than one
            # rotatable bonds.
            if( @allowed_bond_parameters &&
                scalar( @{ $allowed_bond_parameters[0] } ) <
                scalar( keys %bond_parameters ) ) {
                @allowed_bond_parameters =
                    @{ permutation( 2, [],
                                    [ \@allowed_bond_parameters,
                                      \@default_allowed_bond_parameters ], [])};
                @allowed_energies =
                    @{ permutation( 2, [],
                                    [ \@allowed_energies,
                                      \@default_allowed_energies ], [] ) };
                # Flattens angle pairs: [ [ 1 ], [ 2 ] ] =>[ [ 1, 2 ] ].
                @allowed_bond_parameters =
                    map { [ @{ $_->[0] }, @{ $_->[1] } ] }
                        @allowed_bond_parameters;
                @allowed_energies =
                    map { [ $_->[0][0] ] } @allowed_energies;
            } elsif( ! @allowed_bond_parameters ) {
                @allowed_bond_parameters = @default_allowed_bond_parameters;
                @allowed_energies = @default_allowed_energies;
            }

            # Marks visited atoms.
            push @visited_atom_ids, $atom_id;

            # Marks neighbouring atoms.
            push @neighbour_atom_ids,
                @{ $atom_site->{$atom_id}{'connections'} };
            push @neighbour_atom_ids,
                @{ $atom_site->{$atom_id}{'connections_hetatom'} }
            if $include_hetatoms &&
                defined $atom_site->{$atom_id}{'connections_hetatom'};

            # Starts calculating potential energy.
            my ( $next_allowed_bond_parameters, $next_allowed_energies ) =
                @{ threading(
                       \&calc_favourable_angle,
                       { 'parameters' => $parameters,
                         'atom_site' => $atom_site,
                         'atom_id' => $atom_id,
                         'interaction_site' => $interaction_site,
                         'non_bonded_potential' => $non_bonded_potential,
                         'bonded_potential' => $bonded_potential },
                       [ \@allowed_bond_parameters, \@allowed_energies ],
                       $threads ) };

            # my ( $next_allowed_bond_parameters, $next_allowed_energies ) =
            #     @{ calc_favourable_angle(
            #            { 'parameters' => $parameters,
            #              'atom_site' => $atom_site,
            #              'atom_id' => $atom_id,
            #              'interaction_site' => $interaction_site,
            #              'non_bonded_potential' => $non_bonded_potential,
            #              'bonded_potential' => $bonded_potential },
            #            [ \@allowed_bond_parameters, \@allowed_energies ] ) };

            if( scalar @{ $next_allowed_bond_parameters } > 0 ) {
                @allowed_bond_parameters = @{ $next_allowed_bond_parameters };
                @allowed_energies = @{ $next_allowed_energies };
                print info(
                    { message =>
                          $residue_site->{$atom_id}{'pdbx_PDB_model_num'} . " ".
                          $residue_site->{$atom_id}{'label_asym_id'} . " " .
                          $residue_site->{$atom_id}{'label_seq_id'} . " " .
                          $residue_site->{$atom_id}{'label_alt_id'} . " " .
                          "${residue_name} " .
                          $residue_site->{$atom_id}{'label_atom_id'} . " " .
                          "${last_angle_name} " .
                          scalar( @allowed_bond_parameters ) . "\n",
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

    return \@allowed_bond_parameters;
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

    # TODO: not optimal. Angles should be passed.
    my $rotatable_bonds = $atom_site->{$atom_id}{'rotatable_bonds'};
    my $stretchable_bonds = $atom_site->{$atom_id}{'stretchable_bonds'};
    my $bendable_angles = $atom_site->{$atom_id}{'bendable_angles'};

    my %bond_parameters = (
        ( defined $rotatable_bonds ? %{ $rotatable_bonds } : () ),
        ( defined $stretchable_bonds ? %{ $stretchable_bonds } : () ),
        ( defined $bendable_angles ? %{ $bendable_angles } : () )
    );

    my @angle_names =
        sort { $bond_parameters{$a}{'order'} <=> $bond_parameters{$b}{'order'} }
        keys %bond_parameters;

    my @allowed_angles;
    my @allowed_energies;
    for( my $i = 0; $i <= $#{ $array_blocks->[0] }; $i++ ) {
        my $angles = $array_blocks->[0][$i];
        my $energies = $array_blocks->[1][$i][0];
        my %angles =
            map { $angle_names[$_] => [ $angles->[$_] ] }
                ( 0..$#{ $angles } );

        my $pseudo_atom_site =
            generate_pseudo( { 'parameters' => $parameters,
                               'atom_site' => $atom_site,
                               'atom_specifier' => { 'id' => [ "$atom_id" ] },
                               'bond_parameter_values' => \%angles } );
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

    # TODO: not optimal. Angles should be passed.
    my $dihedral_angles =
        collect_dihedral_angles( $residue_site )->{$residue_unique_key};
    my $bendable_angles =
        collect_bond_angles( $residue_site )->{$residue_unique_key};
    my $stretchable_bonds =
        collect_bond_lengths( $residue_site )->{$residue_unique_key};

    my %bond_parameters = (
        ( defined $dihedral_angles ? %{ $dihedral_angles } : () ),
        ( defined $stretchable_bonds ? %{ $stretchable_bonds } : () ),
        ( defined $bendable_angles ? %{ $bendable_angles } : () ),
    );

    my @angle_names =
        sort { $bond_parameters{$a}{'order'} <=> $bond_parameters{$b}{'order'} }
        keys %bond_parameters;

    # Checks for inter-atom interactions and determines if energies
    # comply with cutoffs.
    my @checkable_angles = @{ $array_blocks->[0] };
    my @allowed_angles;
    my @energy_sums;
    my @rmsd_averages;

  ALLOWED_ANGLES:
    for( my $i = 0; $i <= $#checkable_angles; $i++ ) {
        my %rotamer_site = %{ $residue_site };
        my %angles =
            map { ( $angle_names[$_] => $checkable_angles[$i][$_] ) }
                ( 0..$#{ $checkable_angles[$i] } );
        replace_with_rotamer( $parameters, \%rotamer_site, $residue_unique_key,
                              \%angles );

        my @rotamer_atom_ids =
            sort keys %{ filter_new( \%rotamer_site,
                                 { 'exclude' =>
                                   { 'label_atom_id' =>
                                         $interaction_atom_names } } ) };

        # HACK: make sure that $interaction_site atom ids are updated by
        # %rotamer_site.
        my %rotamer_interaction_site = ( %{ $interaction_site }, %rotamer_site );
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

#
# Assigns heteroatoms to specific residues according to "_struct_conn" -- either
# creating the new one or assigning to the existing one. Duplicated heteroatoms
# are to be expected.
# Input:
#     $parameters - general parameters (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $struct_conn - reads 'struc_conn' and assings connections appropriately.
# Output:
#     atom site with assigned heteroatoms.
#

sub assign_hetatoms_with_struct_conn
{
    my ( $parameters, $atom_site, $struct_conn ) = @_;
    $struct_conn //= {};

    return if ! %{ $struct_conn };

    my $hetatom_site =
        filter_new( $atom_site,
                    { 'include' => { 'group_PDB' => [ 'HETATM' ] } } );

    my $last_atom_id = max( keys %{ $atom_site } ) + 1;
    for my $hetatom_id ( keys %{ $hetatom_site } ) {
        my ( $hetatom_label_seq_id, $hetatom_label_asym_id ) =
            map { $hetatom_site->{$hetatom_id}{$_} }
                ( 'label_seq_id',
                  'label_asym_id' );
        my $hetatom_struct_conn =
            filter_new( $struct_conn,
                        { 'include' =>
                          { 'ptnr2_label_seq_id' => [ $hetatom_label_seq_id ],
                            'ptnr2_label_asym_id' => [ $hetatom_label_asym_id ] } } );
        for my $hetatom_struct_conn_id ( sort keys %{ $hetatom_struct_conn } ) {
            my $connected_atom_site =
                filter_new( $atom_site,
                            { 'include' =>
                              { 'label_seq_id' => [
                                    $hetatom_struct_conn->{$hetatom_struct_conn_id}
                                                          {'ptnr1_label_seq_id'}
                                ],
                                'label_asym_id' => [
                                        $hetatom_struct_conn->{$hetatom_struct_conn_id}
                                                              {'ptnr1_label_asym_id'}
                                ],
                                'label_atom_id' => [
                                    $hetatom_struct_conn->{$hetatom_struct_conn_id}
                                                          {'ptnr1_label_atom_id'}
                                ] } } );
            # Iteration has to be performed, because '_struct_conn' does not
            # have 'pdbx_PDB_model_num' and 'label_alt_id' entries.
            my %visited_residues_by_hetatom = ();
            for my $connected_atom_id ( keys %{ $connected_atom_site } ) {
                my $connected_unique_residue_key = unique_residue_key(
                    $connected_atom_site->{$connected_atom_id}
                );

                next if defined $visited_residues_by_hetatom{$connected_unique_residue_key} &&
                    $visited_residues_by_hetatom{$connected_unique_residue_key} == $hetatom_id;
                $visited_residues_by_hetatom{$connected_unique_residue_key} = $hetatom_id;

                # Heteroatom inherits residue information from the atom that is
                # connected to.
                my $current_hetatom = clone $atom_site->{$hetatom_id};
                foreach( 'label_seq_id', 'label_asym_id', 'label_alt_id',
                         'pdbx_PDB_model_num' ) {
                    $current_hetatom->{$_} = $atom_site->{$connected_atom_id}{$_};
                }
                $current_hetatom->{'id'} = $last_atom_id;
                $current_hetatom->{'origin_atom_id'} = $hetatom_id;
                # HACK: Carefully analyse when assigning heteroatom-heteroatom
                # connections.
                if( defined $current_hetatom->{'hybridization'} ) {
                    $current_hetatom->{'hybridization'} = '.';
                }
                $atom_site->{$last_atom_id} = $current_hetatom;

                connect_atoms_explicitly( $atom_site,
                                          [ $last_atom_id ],
                                          [ $connected_atom_id ] );

                $last_atom_id++;
            }
        }

        delete $atom_site->{$hetatom_id};
    }

    return;
}

#
# Assigns heteroatoms to specific residues connecting them to N, O, P, S
# atoms -- either creating the new one or assigning to the existing one.
# Duplicated heteroatoms are to be expected.
# Input:
#     $parameters - general parameters (see Parameters.pm);
#     $atom_site - atom site data structure (see PDBxParser.pm);
# Output:
#     atom site with assigned heteroatoms.
#

sub assign_hetatoms_no_struct_conn
{
    my ( $parameters, $atom_site ) = @_;

    my $hetatom_site =
        filter_new( $atom_site,
                    { 'include' => { 'group_PDB' => [ 'HETATM' ] } } );
    my $interaction_atom_site =
        filter_new( $atom_site,
                    { 'include' =>
                      { 'type_symbol' => [ 'N', 'O', 'P', 'S' ] } } );

    my $last_atom_id = max( keys %{ $atom_site } ) + 1;
    for my $hetatom_id ( sort keys %{ $hetatom_site } ) {
        # NOTE: phosphorus is chosen as max interaction distance as it has the
        # largest vdW radius.
        my $interaction_distance =
            $parameters->{'_[local]_force_field'}{'cutoff_start'} *
            ( ( $parameters->{'_[local]_atom_properties'}
                             {$hetatom_site->{$hetatom_id}{'type_symbol'}}
                             {'vdw_radius'} / 2 ) +
              ( $parameters->{'_[local]_atom_properties'}
                             {'P'}{'vdw_radius'} / 2 ) );
        my $around_site =
            around_distance( $parameters,
                             { $hetatom_id => $hetatom_site->{$hetatom_id},
                               %{ $interaction_atom_site } },
                             { 'id' => [ $hetatom_id ] },
                             $interaction_distance );

        next if ! %{ $around_site };

        for my $around_atom_id ( sort { $a <=> $b } keys %{ $around_site } ) {
            my $around_unique_residue_key = unique_residue_key(
                $around_site->{$around_atom_id}
            );

            # Heteroatom inherits residue information from the atom that is
            # connected to.
            my $current_hetatom = clone $atom_site->{$hetatom_id};
            foreach( 'label_seq_id', 'label_asym_id', 'label_alt_id',
                     'pdbx_PDB_model_num' ) {
                $current_hetatom->{$_} = $atom_site->{$around_atom_id}{$_};
            }
            $current_hetatom->{'id'} = $last_atom_id;
            $current_hetatom->{'origin_atom_id'} = $hetatom_id;
            if( defined $current_hetatom->{'hybridization'} ) {
                $current_hetatom->{'hybridization'} = '.';
            }
            $atom_site->{$last_atom_id} = $current_hetatom;

            connect_atoms_explicitly( $atom_site,
                                      [ $last_atom_id ],
                                      [ $around_atom_id ] );

            $last_atom_id++;
        }

        delete $atom_site->{$hetatom_id};
    }

    return;
}

#
# Quickly determines if the
# Input:
#     $bond_parameters - bond parameter data structure;
#     $bond_parameter_type - bond parameter type;
# Output:
#     if the bond parameter is present or not.
#

sub is_bond_parameter_present
{
    my ( $bond_parameters, $bond_parameter_type ) = @_;
    for my $residue_name ( keys %{ $bond_parameters } ) {
        for my $bond_parameter_name ( keys %{ $bond_parameters->{$residue_name} } ) {
            if( $bond_parameters->{$residue_name}{$bond_parameter_name}{'type'} eq
                $bond_parameter_type ) {
                return 1 if exists $bond_parameters->{$residue_name};
            }
        }
    }
    return 0;
}

1;
