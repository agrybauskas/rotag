package ForceField::Generate;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( obtain_force_field_data
                     force_field_parameters );

use List::Util qw( uniq );

use PDBxParser qw( pdbx_loop_to_array
                   obtain_pdbx_line_new
                   obtain_pdbx_loop );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub obtain_force_field_data
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

sub force_field_parameters
{
    my ( $pdbx_file ) = @_;

    my $force_field_data = obtain_force_field_data( $pdbx_file );

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

        $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                               {'lone_pairs'} = $lone_pair_count;
        $force_field_parameters{'_[local]_atom_properties'}{$type_symbol}
                               {'vdw_radius'} = $vdw_radius;
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

    return;
}

1;
