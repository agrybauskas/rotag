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

    return;
}

1;
