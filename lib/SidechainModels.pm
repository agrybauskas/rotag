package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_translation );

use List::MoreUtils qw( uniq );

use AlterMolecule qw( angle_bending
                      bond_stretching
                      bond_torsion );
use AtomProperties qw( sort_atom_ids_by_name );
use LinearAlgebra qw( mult_matrix_product
                      reshape );
use PDBxParser qw( determine_residue_keys
                   filter_by_unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------- Idealistic sidechain models ------------------------- #

#
# Model that uses rotation around single bonds, angle bending and stretching.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $bond_parameters - bond parameters (see BondParameters.pm);
# Output:
#     adds conformation variable for each atom as list of transformation
#     matrices.
#

sub rotation_translation
{
    my ( $parameters, $atom_site, $bond_parameters, $options ) = @_;

    # Determines all residue ids present in atom site.
    my @residue_unique_keys =
        @{ determine_residue_keys( $atom_site, { 'exclude_dot' => 1 } ) };

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_unique_key ( @residue_unique_keys ) {
        my $residue_site =
            filter_by_unique_residue_key( $atom_site, $residue_unique_key, 1 );

        next if ! %{ $residue_site };

        my $rotatable_bonds = $bond_parameters->rotatable_bonds;
        my $stretchable_bonds = $bond_parameters->stretchable_bonds;
        my $bendable_angles = $bond_parameters->bendable_angles;

        if( ! %{ $rotatable_bonds } &&
            ! %{ $stretchable_bonds } &&
            ! %{ $bendable_angles } ) { next; }

        for my $atom_id ( keys %{ $residue_site }  ) {
            if( ! exists $rotatable_bonds->{$atom_id} &&
                ! exists $stretchable_bonds->{$atom_id} &&
                ! exists $bendable_angles->{$atom_id} ) { next; }

            my @atom_coord = ( $atom_site->{"$atom_id"}{'Cartn_x'},
                               $atom_site->{"$atom_id"}{'Cartn_y'},
                               $atom_site->{"$atom_id"}{'Cartn_z'}, );

            # Matrices for transforming atom coordinates.
            my @bond_torsion_matrices = ();
            my @bond_stretching_matrices = ();
            my @angle_bending_matrices = ();

            if( $do_bond_torsion ) {
                @bond_torsion_matrices =
                    @{ bond_torsion_matrices( $parameters,
                                              $residue_site,
                                              $atom_id,
                                              $rotatable_bonds ) };
            }

            if( $do_bond_stretching ) {
                @bond_stretching_matrices =
                    @{ bond_stretching_matrices( $parameters,
                                                 $residue_site,
                                                 $atom_id,
                                                 $stretchable_bonds ) };
            }

            if( $do_angle_bending ) {
                @angle_bending_matrices =
                    @{ angle_bending_matrices( $parameters,
                                               $residue_site,
                                               $atom_id,
                                               $bendable_angles ) };
            }

            $atom_site->{$atom_id}{'conformation'} =
                mult_matrix_product(
                    [ @bond_torsion_matrices,
                      @angle_bending_matrices,
                      @bond_stretching_matrices,
                      @{ reshape( [ @atom_coord, 1 ], [ 4, 1 ] ) } ] );
        }
    }

    return;
}

sub bond_torsion_matrices
{
    my ( $parameters, $atom_site, $atom_id, $rotatable_bonds ) = @_;

    my @bond_torsion_matrices = ();
    for my $angle_name ( sort { $rotatable_bonds->{$atom_id}{$a}{'order'} <=>
                                $rotatable_bonds->{$atom_id}{$b}{'order'} }
                         keys %{ $rotatable_bonds->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
            map { $rotatable_bonds->{$atom_id}{$angle_name}{'atom_ids'}[$_] }
                ( 2, 1, 0 );
        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $atom_site->{$_}{'Cartn_x'},
                    $atom_site->{$_}{'Cartn_y'},
                    $atom_site->{$_}{'Cartn_z'} ] }
                ( $mid_atom_id, $up_atom_id, $side_atom_id );

        push @bond_torsion_matrices,
            @{ bond_torsion( $parameters,
                             $mid_atom_coord,
                             $up_atom_coord,
                             $side_atom_coord,
                             $angle_name ) };
    }

    return \@bond_torsion_matrices;
}

sub bond_stretching_matrices
{
    my ( $parameters, $atom_site, $atom_id, $stretchable_bonds ) = @_;

    # TODO: it is possible to cache here.
    my @bond_stretching_matrices = ();
    for my $bond_name ( sort { $stretchable_bonds->{$atom_id}{$a}{'order'} <=>
                               $stretchable_bonds->{$atom_id}{$b}{'order'} }
                        keys %{ $stretchable_bonds->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id ) =
            map { $stretchable_bonds->{$atom_id}{$bond_name}{'atom_ids'}[$_] }
                ( 1, 0 );
        my @mid_connections = # Excludes up atom.
            grep { $_ ne $up_atom_id }
                @{ $atom_site->{$mid_atom_id}{'connections'} };
        my ( $side_atom_id ) =
            @{ sort_atom_ids_by_name( \@mid_connections, $atom_site ) };

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $atom_site->{$_}{'Cartn_x'},
                    $atom_site->{$_}{'Cartn_y'},
                    $atom_site->{$_}{'Cartn_z'} ] }
                ( $mid_atom_id, $up_atom_id, $side_atom_id );

        push @bond_stretching_matrices,
             @{ bond_stretching( $parameters,
                                 $mid_atom_coord,
                                 $up_atom_coord,
                                 $side_atom_coord,
                                 $bond_name ) };
    }

    return \@bond_stretching_matrices;
}

sub angle_bending_matrices
{
    my ( $parameters, $atom_site, $atom_id, $bendable_angles ) = @_;

    my @angle_bending_matrices = ();
    for my $angle_name ( sort { $bendable_angles->{$atom_id}{$a}{'order'} <=>
                                $bendable_angles->{$atom_id}{$b}{'order'} }
                         keys %{ $bendable_angles->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
            map { $bendable_angles->{$atom_id}{$angle_name}{'atom_ids'}[$_] }
                ( 2, 1, 0 );

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $atom_site->{$_}{'Cartn_x'},
                    $atom_site->{$_}{'Cartn_y'},
                    $atom_site->{$_}{'Cartn_z'} ] }
                ( $mid_atom_id, $up_atom_id, $side_atom_id );

        push @angle_bending_matrices,
             @{ angle_bending( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               $angle_name,
                               'eta' ) }; # 'eta' should be equal to 0.
    }

    return \@angle_bending_matrices;
}

1;
