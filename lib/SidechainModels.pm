package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_only
                     rotation_translation );

use List::MoreUtils qw( uniq );

use AlterMolecule qw( angle_bending
                      bond_stretching
                      bond_torsion );
use AtomProperties qw( sort_atom_names );
use BondProperties qw( bendable_angles
                       rotatable_bonds
                       stretchable_bonds );
use LinearAlgebra qw( mult_matrix_product
                      reshape );
use PDBxParser qw( determine_residue_keys
                   filter_new
                   filter_by_unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------- Idealistic sidechain models ------------------------- #

#
# Model that uses only rotation around single bonds.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure.
# Output:
#     adds conformation variable for each atom as list of transformation
#     matrices.
#

sub rotation_only
{
    my ( $parameters, $atom_site ) = @_;
    return rotation_translation( $parameters, $atom_site,
                                 { 'do_bond_torsion' => 1,
                                   'do_bond_stretching' => 0,
                                   'do_angle_bending' => 0 } );
}

#
# Model that uses rotation around single bonds, angle bending and stretching.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $options{'do_bond_torsion'} - calculates bond torsion matrices;
#     $options{'do_bond_stretching'} - calculates bond stretching  matrices;
#     $options{'do_angle_bending'} - calculates angle bending matrices;
#     $options{'include_hetatoms'} - includes heteroatoms to the calculations.
#
# Output:
#     adds conformation variable for each atom as list of transformation
#     matrices.
#

sub rotation_translation
{
    my ( $parameters, $atom_site, $options ) = @_;
    my ( $do_bond_torsion, $do_bond_stretching, $do_angle_bending,
         $include_hetatoms ) = (
        $options->{do_bond_torsion},
        $options->{do_bond_stretching},
        $options->{do_angle_bending},
        $options->{include_hetatoms},
    );

    $do_bond_torsion //= 1;
    $do_bond_stretching //= 1;
    $do_angle_bending //= 1;
    $include_hetatoms //= 0;

    my %atom_site = %{ $atom_site }; # Copy of $atom_site.

    # Determines all residue ids present in atom site.
    my @residue_unique_keys =
        @{ determine_residue_keys( $atom_site, { 'exclude_dot' => 1 } ) };

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_unique_key ( @residue_unique_keys ) {
        my $residue_site =
            filter_by_unique_residue_key( \%atom_site, $residue_unique_key, 1 );

        next if ! %{ $residue_site };

        my $bendable_angles =
            $do_angle_bending ?
            bendable_angles( $parameters, $residue_site, undef,
                             { 'include_hetatoms' => $include_hetatoms } ) : {};
        my $rotatable_bonds =
            $do_bond_torsion ?
            rotatable_bonds( $parameters, $residue_site, undef,
                             { 'include_hetatoms' => $include_hetatoms } ) : {};
        my $stretchable_bonds =
            $do_bond_stretching ?
            stretchable_bonds( $parameters, $residue_site, undef,
                               { 'include_hetatoms' => $include_hetatoms } ): {};

        if( ! %{ $rotatable_bonds } &&
            ! %{ $stretchable_bonds } &&
            ! %{ $bendable_angles } ) { next; }

        for my $atom_id ( keys %{ $residue_site }  ) {
            if( ! exists $rotatable_bonds->{$atom_id} &&
                ! exists $stretchable_bonds->{$atom_id} &&
                ! exists $bendable_angles->{$atom_id} ) { next; }

            my @atom_coord = ( $atom_site{"$atom_id"}{'Cartn_x'},
                               $atom_site{"$atom_id"}{'Cartn_y'},
                               $atom_site{"$atom_id"}{'Cartn_z'}, );

            # Matrices for transforming atom coordinates.
            my @bond_torsion_matrices;
            my @bond_stretching_matrices;
            my @angle_bending_matrices;

            if( $do_bond_torsion ) {
                # TODO: code block is similar to all_dihedral(). The code should
                # be moved to separate function.
                for my $angle_name (
                    sort { $rotatable_bonds->{$atom_id}{$a}{'order'} <=>
                           $rotatable_bonds->{$atom_id}{$b}{'order'} }
                    keys %{ $rotatable_bonds->{$atom_id} } ) {
                    my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
                        map { $rotatable_bonds->{$atom_id}{$angle_name}{'atom_ids'}[$_] }
                            ( 2, 1, 0 );

                    my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
                        map { [ $residue_site->{$_}{'Cartn_x'},
                                $residue_site->{$_}{'Cartn_y'},
                                $residue_site->{$_}{'Cartn_z'} ] }
                            ( $mid_atom_id, $up_atom_id, $side_atom_id );

                    # Creates and appends matrices to a list of matrices that
                    # later will be multiplied.
                    push @bond_torsion_matrices,
                         @{ bond_torsion( $parameters,
                                          $mid_atom_coord,
                                          $up_atom_coord,
                                          $side_atom_coord,
                                          $angle_name ) };
                }
            }

            if( $do_bond_stretching ) {
                for my $bond_name (
                    sort { $stretchable_bonds->{$atom_id}{$a}{'order'} <=>
                           $stretchable_bonds->{$atom_id}{$b}{'order'} }
                    keys %{ $stretchable_bonds->{$atom_id} } ) {
                    my $up_atom_id =
                        $stretchable_bonds->{$atom_id}{$bond_name}{'atom_ids'}[1];
                    my $mid_atom_id =
                        $stretchable_bonds->{$atom_id}{$bond_name}{'atom_ids'}[0];

                    my @mid_connections = # Excludes up atom.
                        grep { $_ ne $up_atom_id }
                            @{ $residue_site->{$mid_atom_id}{'connections'} };
                    my @mid_connection_names = # Excludes up atom.
                        map { $residue_site->{$_}{'label_atom_id'} }
                            @mid_connections;
                    my $side_atom_name =
                        sort_atom_names( \@mid_connection_names )->[0];
                    my $side_atom_id =
                        filter_new( $residue_site,
                                {  'include' =>
                                       { 'label_atom_id' => [ $side_atom_name ] },
                                   'return_data' => 'id' } )->[0];

                    my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
                        map { [ $residue_site->{$_}{'Cartn_x'},
                                $residue_site->{$_}{'Cartn_y'},
                                $residue_site->{$_}{'Cartn_z'} ] }
                            ( $mid_atom_id, $up_atom_id, $side_atom_id );

                    # Creates and appends matrices to a list of matrices that
                    # later will be multiplied.
                    push @bond_stretching_matrices,
                         @{ bond_stretching( $parameters,
                                             $mid_atom_coord,
                                             $up_atom_coord,
                                             $side_atom_coord,
                                             $bond_name ) };
                }
            }

            if( $do_angle_bending ) {
                for my $angle_name ( sort { $bendable_angles->{$atom_id}{$a}{'order'} <=>
                                            $bendable_angles->{$atom_id}{$b}{'order'} }
                                     keys %{ $bendable_angles->{$atom_id} } ) {
                    my $up_atom_id =
                        $bendable_angles->{$atom_id}{$angle_name}{'atom_ids'}[2];
                    my $mid_atom_id =
                        $bendable_angles->{$atom_id}{$angle_name}{'atom_ids'}[1];
                    my $side_atom_id =
                        $bendable_angles->{$atom_id}{$angle_name}{'atom_ids'}[0];

                    my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
                        map { [ $residue_site->{$_}{'Cartn_x'},
                                $residue_site->{$_}{'Cartn_y'},
                                $residue_site->{$_}{'Cartn_z'} ] }
                            ( $mid_atom_id, $up_atom_id, $side_atom_id );

                    # Creates and appends matrices to a list of matrices that
                    # later will be multiplied.
                    push @angle_bending_matrices,
                         @{ angle_bending( $parameters,
                                           $mid_atom_coord,
                                           $up_atom_coord,
                                           $side_atom_coord,
                                           $angle_name,
                                           'eta' ) }; # 'eta' should equal to 0.
                }
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

1;
