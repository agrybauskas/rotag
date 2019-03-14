package Parameters;

use strict;
use warnings;

use File::Basename qw( dirname );
use List::Util qw( max
                   uniq );

use PDBxParser qw( pdbx_loop_to_array
                   obtain_pdbx_line_new
                   obtain_pdbx_loop );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my ( $force_field_file ) = ( $args->{'force_field_file'} );

    $force_field_file //= dirname( __FILE__ ) . '/Parameters.cif';

    my $force_field = force_field( $force_field_file );
    my $constants = constants( $force_field );
    my $bond_types = bond_types( $force_field );
    my $self = { %{ $constants }, %{ $bond_types }, %{ $force_field } };

    return bless $self, $class;
}

# ----------------------------- Member functions ------------------------------ #

sub _obtain_force_field_data
{
    my ( $pdbx_file ) = @_;

    my @looped_categories;
    my @unlooped_items;

    my $is_looped = 0;

    local @ARGV = ( $pdbx_file );
    while( <> ) {
        if( m/^loop_/ ) {
            $is_looped = 1;
        } elsif( $is_looped && m/^(_\S+)\./ ) {
            my $category = $1;
            push @looped_categories, $category;
        } elsif( m/(^_\S+\.\S+)/ ) {
            my $item = $1;
            push @unlooped_items, $item;
        } else {
            $is_looped = 0;
        }
    }

    @looped_categories = uniq @looped_categories;

    my $pdbx_line_data = obtain_pdbx_line_new( $pdbx_file, \@unlooped_items );
    my $pdbx_loop_data = obtain_pdbx_loop( $pdbx_file, \@looped_categories );

    return { %{ $pdbx_line_data }, %{ $pdbx_loop_data } };
}

# sub _covalent_bond_combinations
# {
#     my ( $ATOMS ) = @_;

#     my %covalent_bond_combinations;
#     for my $atom_i_name ( keys %Parameters::ATOMS ) {
#         for my $atom_j_name ( keys %Parameters::ATOMS ) {
#             $covalent_bond_combinations{$atom_i_name}{$atom_j_name}{'length'} =
#                 permutation( 2,
#                              [],
#                              [ $Parameters::ATOMS{$atom_i_name}{'covalent_radius'}{'length'},
#                                $Parameters::ATOMS{$atom_j_name}{'covalent_radius'}{'length'} ],
#                              [] );
#             $covalent_bond_combinations{$atom_i_name}{$atom_j_name}{'error'} =
#                 permutation( 2,
#                              [],
#                              [ $Parameters::ATOMS{$atom_i_name}{'covalent_radius'}{'error'},
#                                $Parameters::ATOMS{$atom_j_name}{'covalent_radius'}{'error'} ],
#                              [] );
#         }
#     }

#     return \%covalent_bond_combinations;
# }

sub _epsilon
{
    my $EPSILON = 1.0;

    while( ( 1.0 + 0.5 * $EPSILON ) != 1.0 ) {
        $EPSILON = 0.5 * $EPSILON;
    }

    return $EPSILON;
}

sub _pi
{
    return 4 * atan2 1, 1;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub constants
{
    my ( $force_field ) = @_;

    my %constants;

    # General constants.
    $constants{'_constants'}{'epsilon'} = _epsilon();
    $constants{'_constants'}{'pi'} = _pi();
    $constants{'_constants'}{'sig_figs_min'} = '%.3f';
    $constants{'_constants'}{'sig_figs_max'} = '%.6f';

    # Angle constants.
    $constants{'_constants'}{'sp3_angle'} =
        109.5 * $constants{'_constants'}{'pi'} / 180.0;
    $constants{'_constants'}{'sp2_angle'} =
        120.0 * $constants{'_constants'}{'pi'} / 180.0;
    $constants{'_constants'}{'sp2_angle'} =
        $constants{'_constants'}{'pi'};
    $constants{'_constants'}{'sp3_angle_err'} =
        5.0 * $constants{'_constants'}{'pi'} / 180.0;
    $constants{'_constants'}{'sp2_angle_err'} =
        5.0 * $constants{'_constants'}{'pi'} / 180.0;
    $constants{'_constants'}{'sp_angle_err'} =
        5.0 * $constants{'_constants'}{'pi'} / 180.0;

    # Grid constants.
    $constants{'_constants'}{'edge_length_connection'} =
        max( map { @{ $force_field->{'_[local]_atom_properties'}{$_}
                                    {'covalent_radius'}{'length'} } }
             keys %{ $force_field->{'_[local]_atom_properties'} } ) * 2;

    # Arginine model is use for calculating the interaction cutoff.
    $constants{'_constants'}{'edge_length_interaction'} =
        7 * $force_field->{'_[local]_atom_properties'}{'C'}{'covalent_radius'}
                          {'length'}->[0] +
        2 * $force_field->{'_[local]_atom_properties'}{'N'}{'covalent_radius'}
                          {'length'}->[0] +
        3 * $force_field->{'_[local]_atom_properties'}{'C'}{'covalent_radius'}
                          {'length'}->[1] +
        3 * $force_field->{'_[local]_atom_properties'}{'N'}{'covalent_radius'}
                          {'length'}->[1] +
            $force_field->{'_[local]_atom_properties'}{'H'}{'covalent_radius'}
                          {'length'}->[0];

    return \%constants;
}

sub force_field
{
    my ( $pdbx_file ) = @_;

    my $force_field_data = _obtain_force_field_data( $pdbx_file );

    my %force_field_parameters = ();

    # Restructuring force field constants.
    $force_field_parameters{'_[local]_force_field'} =
        pdbx_loop_to_array( $force_field_data, '_[local]_force_field' );

    # Restructuring parameters of atom properties.
    my $atom_properties_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_atom_properties' );

    for my $atom_properties ( @{ $atom_properties_loop } ) {
        my $type_symbol = $atom_properties->{'type_symbol'};
        my $lone_pair_count = $atom_properties->{'lone_pair_count'};
        my $vdw_radius = $atom_properties->{'vdw_radius'};
        my $valence = $atom_properties->{'valence'};
        my $covalent_radius_value = $atom_properties->{'covalent_radius_value'};
        my $covalent_radius_error = $atom_properties->{'covalent_radius_error'};

        $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                               {'lone_pairs'} = $lone_pair_count;
        $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                               {'vdw_radius'} = $vdw_radius;
        $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                               {'valence'} = $valence;

        push @{ $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                                       {'covalent_radius'}{'length'} },
            $covalent_radius_value;
        push @{ $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                                       {'covalent_radius'}{'error'} },
            $covalent_radius_error;
    }

    # Restructuring parameters of Lennard-Jones.
    my $lennard_jones_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_lennard_jones' );

    for my $lennard_jones ( @{ $lennard_jones_loop } ) {
        my $type_symbol_1 = $lennard_jones->{'type_symbol_1'};
        my $type_symbol_2 = $lennard_jones->{'type_symbol_2'};
        my $sigma = $lennard_jones->{'sigma'};
        my $epsilon = $lennard_jones->{'epsilon'};

        $force_field_parameters{'_[local]_lennard_jones'}{$type_symbol_1}
                               {$type_symbol_2}{'sigma'} = $sigma;
        $force_field_parameters{'_[local]_lennard_jones'}{$type_symbol_1}
                               {$type_symbol_2}{'epsilon'} = $epsilon;
    }

    # Restructuring parameters of partial charge.
    my $partial_charge_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_partial_charge' );

    for my $partial_charge ( @{ $partial_charge_loop } ) {
        my $residue_name = $partial_charge->{'label_comp_id'};
        my $atom_name = $partial_charge->{'label_atom_id'};
        my $partial_charge_value = $partial_charge->{'value'};

        $force_field_parameters{'_[local]_partial_charge'}{$residue_name}
                               {$atom_name} = $partial_charge_value;
    }

    # Restructuring parameters of partial torsional potential.
    my $torsional_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_torsional' );

    for my $torsional ( @{ $torsional_loop } ) {
        my $type_symbol_1 = $torsional->{'type_symbol_1'};
        my $type_symbol_2 = $torsional->{'type_symbol_2'};
        my $epsilon = $torsional->{'epsilon'};

        $force_field_parameters{'_[local]_torsional'}{$type_symbol_1}
                               {$type_symbol_2}{'epsilon'} = $epsilon;
    }

    # Restructuring parameters of hydrogen bond.
    my $hydrogen_bond_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_h_bond' );

    for my $hydrogen_bond ( @{ $hydrogen_bond_loop } ) {
        my $type_symbol = $hydrogen_bond->{'type_symbol'};
        my $sigma = $hydrogen_bond->{'sigma'};
        my $epsilon = $hydrogen_bond->{'epsilon'};

        $force_field_parameters{'_[local]_h_bond'}{$type_symbol}{'sigma'} =
            $sigma;
        $force_field_parameters{'_[local]_h_bond'}{$type_symbol}{'epsilon'} =
            $epsilon;
    }

    # Restructuring parameters of atom necessity.
    my $residue_atoms_loop =
        pdbx_loop_to_array($force_field_data, '_[local]_residue_atom_necessity');

    for my $residue_atoms ( @{ $residue_atoms_loop } ) {
        my $residue_name = $residue_atoms->{'label_comp_id'};
        my $atom_name = $residue_atoms->{'label_atom_id'};
        my $necessity_value = $residue_atoms->{'value'};

        push @{ $force_field_parameters{'_[local]_residue_atom_necessity'}
                                       {$residue_name}
                                       {$necessity_value} }, $atom_name;
    }

    # Restructuring parameters of clear hybridizations.
    my $clear_hybridization_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_clear_hybridization' );

    for my $clear_hybridization ( @{ $clear_hybridization_loop } ) {
        my $residue_name = $clear_hybridization->{'label_comp_id'};
        my $atom_name = $clear_hybridization->{'label_atom_id'};
        my $hybridization = $clear_hybridization->{'type'};

        $force_field_parameters{'_[local]_clear_hybridization'}{$residue_name}
                               {$atom_name} = $hybridization;
    }

    # Restructuring parameters of connectivity.
    my $connectivity_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_connectivity' );

    for my $connectivity ( @{ $connectivity_loop } ) {
        my $residue_name = $connectivity->{'label_comp_id'};
        my $atom_name_1 = $connectivity->{'label_atom_1_id'};
        my $atom_name_2 = $connectivity->{'label_atom_2_id'};

        push @{ $force_field_parameters{'_[local]_hydrogen_names'}{$residue_name}
                                       {$atom_name_1} }, $atom_name_2;
    }

    # Restructuring parameters of hydrogen names.
    my $hydrogen_names_loop =
        pdbx_loop_to_array( $force_field_data, '_[local]_hydrogen_names' );

    for my $hydrogen_names ( @{ $hydrogen_names_loop } ) {
        my $residue_name = $hydrogen_names->{'label_comp_id'};
        my $atom_name = $hydrogen_names->{'label_atom_id'};
        my $hydrogen_name = $hydrogen_names->{'label_hydrogen_atom_id'};

        $force_field_parameters{'_[local]_hydrogen_names'}{$residue_name}
                               {$atom_name} = $hydrogen_name;
    }

    # Restructuring parameters of interaction atom names.
    $force_field_parameters{'_[local]_interaction_atoms'} =
        $force_field_data->{'_[local]_interaction_atoms'}{'data'};

    # Restructuring parameters of mainchain atom names.
    $force_field_parameters{'_[local]_mainchain_atom_names'} =
        $force_field_data->{'_[local]_mainchain_atom_names'}{'data'};

    # Restructuring parameters of sidechain atom names.
    $force_field_parameters{'_[local]_sidechain_atom_names'} =
        $force_field_data->{'_[local]_sidechain_atom_names'}{'data'};

    # Restructuring parameters of rotatable residue names.
    $force_field_parameters{'_[local]_rotatable_residue_names'} =
        $force_field_data->{'_[local]_rotatable_residue_names'}{'data'};

    return \%force_field_parameters;
}

sub bond_types
{
    my ( $force_field ) = @_;

    my %bond_types;

    # Single bonds.
    my @atom_symbols_single = qw( H C N O S );
    for my $first_atom_symbol ( @atom_symbols_single ) {
        for my $second_atom_symbol ( @atom_symbols_single ) {
            $bond_types{'_bond_types'}{'single'}{$first_atom_symbol}
                       {$second_atom_symbol}{'min_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[0] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[0] -
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[0] -
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[0];
            $bond_types{'_bond_types'}{'single'}{$first_atom_symbol}
                       {$second_atom_symbol}{'max_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[0] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[0] +
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[0] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[0];
        }
    }

    # Double bonds.
    my @atom_symbols_double = qw( C N O );
    for my $first_atom_symbol ( @atom_symbols_double ) {
        for my $second_atom_symbol ( @atom_symbols_double ) {
            $bond_types{'_bond_types'}{'double'}{$first_atom_symbol}
                       {$second_atom_symbol}{'min_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[1] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[1] -
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[1] -
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[1];
            $bond_types{'_bond_types'}{'double'}{$first_atom_symbol}
                       {$second_atom_symbol}{'max_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[1] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[1] +
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[1] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[1];
        }
    }

    # Triple bonds.
    my @atom_symbols_triple = qw( C );
    for my $first_atom_symbol ( @atom_symbols_triple ) {
        for my $second_atom_symbol ( @atom_symbols_triple ) {
            $bond_types{'_bond_types'}{'triple'}{$first_atom_symbol}
                       {$second_atom_symbol}{'min_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[2] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[2] -
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[2] -
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[2];
            $bond_types{'_bond_types'}{'triple'}{$first_atom_symbol}
                       {$second_atom_symbol}{'max_length'} =
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'length'}[2] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'length'}[2] +
                $force_field->{'_[local]_atom_properties'}{$first_atom_symbol}
                              {'covalent_radius'}{'error'}[2] +
                $force_field->{'_[local]_atom_properties'}{$second_atom_symbol}
                              {'covalent_radius'}{'error'}[2];
        }
    }

    return \%bond_types;
}

# our %COVALENT_BOND_COMB = %{ covalent_bond_combinations( \%Parameters::ATOMS ) };


1;
