package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_translation
                     rotation_translation_new );

use List::Util qw( uniq );

use AlterMolecule qw( angle_bending
                      bond_altering
                      bond_stretching
                      bond_torsion );
use AtomProperties qw( sort_atom_ids_by_name );
use BondParameters qw( collect_bond_angles
                       collect_bond_lengths
                       collect_dihedral_angles );
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
    my ( $parameters, $atom_site, $options ) = @_;

    $options //= { 'ref_atom_site' => $atom_site };

    # Determines all residue ids present in atom site.
    my @residue_unique_keys =
        @{ determine_residue_keys( $atom_site, { 'exclude_dot' => 1 } ) };

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_unique_key ( @residue_unique_keys ) {
        my $residue_site = filter_by_unique_residue_key(
            $atom_site, $residue_unique_key, 1,
        );

        next if ! %{ $residue_site };

        my $rotatable_bonds =
            { map { $_ => $residue_site->{$_}{'rotatable_bonds'} }
              grep { defined $residue_site->{$_}{'rotatable_bonds'} }
              keys %{ $residue_site } };
        my $stretchable_bonds =
            { map { $_ => $residue_site->{$_}{'stretchable_bonds'} }
              grep { defined $residue_site->{$_}{'stretchable_bonds'} }
              keys %{ $residue_site } };
        my $bendable_angles =
            { map { $_ => $residue_site->{$_}{'bendable_angles'} }
              grep { defined $residue_site->{$_}{'bendable_angles'} }
              keys %{ $residue_site } };

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

            if( %{ $rotatable_bonds } ) {
                @bond_torsion_matrices =
                    @{ bond_torsion_matrices( $parameters,
                                              $residue_site,
                                              $atom_id,
                                              $rotatable_bonds,
                                              $options ) };
            }

            if( %{ $stretchable_bonds } ) {
                @bond_stretching_matrices =
                    @{ bond_stretching_matrices( $parameters,
                                                 $residue_site,
                                                 $atom_id,
                                                 $stretchable_bonds,
                                                 $options ) };
            }

            if( %{ $bendable_angles } ) {
                @angle_bending_matrices =
                    @{ angle_bending_matrices( $parameters,
                                               $residue_site,
                                               $atom_id,
                                               $bendable_angles,
                                               $options ) };
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

sub rotation_translation_new
{
    my ( $parameters, $atom_site ) = @_;

    # Determines all residue ids present in atom site.
    my @residue_unique_keys =
        @{ determine_residue_keys( $atom_site, { 'exclude_dot' => 1 } ) };

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_unique_key ( @residue_unique_keys ) {
        my $residue_site = filter_by_unique_residue_key(
            $atom_site, $residue_unique_key, 1,
        );

        next if ! %{ $residue_site };

        my $rotatable_bonds =
            { map { $_ => $residue_site->{$_}{'rotatable_bonds'} }
              grep { defined $residue_site->{$_}{'rotatable_bonds'} }
              keys %{ $residue_site } };
        my $stretchable_bonds =
            { map { $_ => $residue_site->{$_}{'stretchable_bonds'} }
              grep { defined $residue_site->{$_}{'stretchable_bonds'} }
              keys %{ $residue_site } };
        my $bendable_angles =
            { map { $_ => $residue_site->{$_}{'bendable_angles'} }
              grep { defined $residue_site->{$_}{'bendable_angles'} }
              keys %{ $residue_site } };

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
            my @conformation_matrices =
                @{ conformation_matrices( $parameters,
                                          $residue_site,
                                          $atom_id,
                                          $stretchable_bonds,
                                          $bendable_angles,
                                          $rotatable_bonds,
                                          { 'ref_atom_site' => $atom_site } ) };

            $atom_site->{$atom_id}{'conformation'} =
                mult_matrix_product(
                    [ @conformation_matrices,
                      @{ reshape( [ @atom_coord, 1 ], [ 4, 1 ] ) } ] );
        }
    }

    return;
}

#
# Creates bond torsion matrices that are used to transform atom coordinates.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $atom_id - atom id;
#     $rotatable_bonds - rotatable bonds data structure in BondParameters.pm;
# Output:
#     \@bond_torsion_matrices - bond torsion matrices.
#

sub bond_torsion_matrices
{
    my ( $parameters, $atom_site, $atom_id, $rotatable_bonds, $options ) = @_;
    my ( $ref_atom_site ) = ( $options->{'ref_atom_site'} );
    $ref_atom_site //= $atom_site;

    my @bond_torsion_matrices = ();
    for my $angle_name ( sort { $rotatable_bonds->{$atom_id}{$a}{'order'} <=>
                                $rotatable_bonds->{$atom_id}{$b}{'order'} }
                         keys %{ $rotatable_bonds->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
            map { $rotatable_bonds->{$atom_id}{$angle_name}{'atom_ids'}[$_] }
                ( 2, 1, 0 );

        next if ! exists $ref_atom_site->{$mid_atom_id} ||
                ! exists $ref_atom_site->{$up_atom_id} ||
                ! exists $ref_atom_site->{$side_atom_id};

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $ref_atom_site->{$_}{'Cartn_x'},
                    $ref_atom_site->{$_}{'Cartn_y'},
                    $ref_atom_site->{$_}{'Cartn_z'} ] }
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

#
# Creates bond stretching matrices that are used to transform atom coordinates.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $atom_id - atom id;
#     $stretchable_bonds - stretchable bonds data structure in BondParameters.pm.
# Output:
#     \@bond_stretching_matrices - bond stretching matrices.
#

sub bond_stretching_matrices
{
    my ( $parameters, $atom_site, $atom_id, $stretchable_bonds, $options ) = @_;
    my ( $ref_atom_site ) = ( $options->{'ref_atom_site'} );
    $ref_atom_site //= $atom_site;

    # TODO: it is possible to cache here.
    my @bond_stretching_matrices = ();
    for my $bond_name ( sort { $stretchable_bonds->{$atom_id}{$a}{'order'} <=>
                               $stretchable_bonds->{$atom_id}{$b}{'order'} }
                        keys %{ $stretchable_bonds->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id ) =
            map { $stretchable_bonds->{$atom_id}{$bond_name}{'atom_ids'}[$_] }
                ( 1, 0 );

        next if ! exists $ref_atom_site->{$mid_atom_id} ||
                ! exists $ref_atom_site->{$up_atom_id};

        my @mid_connections = # Excludes up atom.
            grep { exists $ref_atom_site->{$_} }
            grep { $_ ne $up_atom_id }
                 ( @{ $ref_atom_site->{$mid_atom_id}{'connections'} },
                   @{ $ref_atom_site->{$mid_atom_id}{'connections_hetatom'} } );

        my ( $side_atom_id ) =
            @{ sort_atom_ids_by_name( \@mid_connections, $ref_atom_site ) };

        next if ! defined $side_atom_id ||
                ! exists $ref_atom_site->{$side_atom_id};

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $ref_atom_site->{$_}{'Cartn_x'},
                    $ref_atom_site->{$_}{'Cartn_y'},
                    $ref_atom_site->{$_}{'Cartn_z'} ] }
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

#
# Creates angle bending matrices that are used to transform atom coordinates.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $atom_id - atom id;
#     $bendable_angles - bendable angles data structure in BondParameters.pm.
# Output:
#     \@angle_bending_matrices - bond stretching matrices.
#

sub angle_bending_matrices
{
    my ( $parameters, $atom_site, $atom_id, $bendable_angles, $options ) = @_;
    my ( $ref_atom_site ) = ( $options->{'ref_atom_site'} );
    $ref_atom_site //= $atom_site;

    my @angle_bending_matrices = ();
    for my $angle_name ( sort { $bendable_angles->{$atom_id}{$a}{'order'} <=>
                                $bendable_angles->{$atom_id}{$b}{'order'} }
                         keys %{ $bendable_angles->{$atom_id} } ) {

        my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
            map { $bendable_angles->{$atom_id}{$angle_name}{'atom_ids'}[$_] }
                ( 2, 1, 0 );

        next if ! exists $ref_atom_site->{$mid_atom_id} ||
                ! exists $ref_atom_site->{$up_atom_id} ||
                ! exists $ref_atom_site->{$side_atom_id};

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $ref_atom_site->{$_}{'Cartn_x'},
                    $ref_atom_site->{$_}{'Cartn_y'},
                    $ref_atom_site->{$_}{'Cartn_z'} ] }
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

#
# Creates generalized series of matrices that are used to transform atom
# coordinates by changing bond length, bond and dihedral angles.
# Input:
#     $parameters - force-field parameters (see Parameters.pm);
#     $atom_site - atom data structure;
#     $atom_id - atom id;
#     $stretchable_bonds - stretchable bonds data structure in BondParameters.pm;
#     $bendable_angles - bendable angles data structure in BondParameters.pm;
#     $rotatable_bonds - rotatable bonds data structure in BondParameters.pm.
# Output:
#     \@conformation_matrices - bond stretching matrices.
#

sub conformation_matrices
{
    my ( $parameters, $atom_site, $atom_id, $stretchable_bonds, $bendable_angles,
         $rotatable_bonds, $options ) = @_;
    my ( $ref_atom_site ) = ( $options->{'ref_atom_site'} );
    $ref_atom_site //= $atom_site;

    # Sorting by the priority: bond length > bond angle > dihedral angle, because
    # the other parameters can be derived from the previous parameter.
    my ( $bond_parameters ) =
        grep { %{ $_ } }
             ( $stretchable_bonds, $bendable_angles, $rotatable_bonds );

    my @conformation_matrices = ();
    for my $parameter_name ( sort { $bond_parameters->{$atom_id}{$a}{'order'} <=>
                                    $bond_parameters->{$atom_id}{$b}{'order'} ||
                                    $bond_parameters->{$atom_id}{$b}{'rank'} <=>
                                    $bond_parameters->{$atom_id}{$a}{'rank'} }
                             keys %{ $bond_parameters->{$atom_id} } ) {
        my $parameter_data = $bond_parameters->{$atom_id}{$parameter_name};
        my $atom_ids = $parameter_data->{'atom_ids'};

        my @ref_frame_atom_idxs =
            scalar @{ $atom_ids } > 2 ? ( 2, 1, 0 ) : ( 1, 0 );
        my ( $up_atom_id, $mid_atom_id, $side_atom_id ) =
            map { $atom_ids->[$_] } @ref_frame_atom_idxs;

        if( ! defined $side_atom_id ) {
            my @mid_connections =
                uniq
                grep { $_ ne $up_atom_id }
                     ( @{ $ref_atom_site->{$mid_atom_id}{'connections'} } );
            my @mid_heteroatom_connections =
                uniq
                grep { $_ ne $up_atom_id }
                     ( @{ $ref_atom_site->{$mid_atom_id}{'connections_hetatom'} } );
            if( @mid_connections ) {
                ( $side_atom_id ) =
                    @{ sort_atom_ids_by_name( \@mid_connections, $ref_atom_site ) };
            } else {
                ( $side_atom_id ) =
                    @{ sort_atom_ids_by_name( \@mid_heteroatom_connections, $ref_atom_site ) };
            }
        }

        my $bond_name;
        my $bond_angle_name;
        my $dihedral_angle_name;

        next if ! defined $mid_atom_id ||
                ! defined $up_atom_id ||
                ! defined $side_atom_id;

        if( defined $stretchable_bonds && exists $stretchable_bonds->{$atom_id} ) {
            ( $bond_name  ) =
                grep { join( ',', @{ $stretchable_bonds->{$atom_id}{$_}{'atom_ids'} } )
                           =~ m/^${mid_atom_id},${up_atom_id}$/ }
                keys %{ $stretchable_bonds->{$atom_id} };
        }

        if( defined $bendable_angles && exists $bendable_angles->{$atom_id} ) {
            ( $bond_angle_name ) =
                grep { join( ',', @{ $bendable_angles->{$atom_id}{$_}{'atom_ids'} } )
                           =~ m/^${side_atom_id},${mid_atom_id},${up_atom_id}$/ }
                keys %{ $bendable_angles->{$atom_id} };
        }

        if( defined $rotatable_bonds && exists $rotatable_bonds->{$atom_id} ) {
            ( $dihedral_angle_name ) =
                grep { join( ',', @{ $rotatable_bonds->{$atom_id}{$_}{'atom_ids'} } )
                           =~ m/^\d+,${mid_atom_id},${up_atom_id},\d+$/ }
                keys %{ $rotatable_bonds->{$atom_id} };
        }

        my ( $mid_atom_coord, $up_atom_coord, $side_atom_coord ) =
            map { [ $ref_atom_site->{$_}{'Cartn_x'},
                    $ref_atom_site->{$_}{'Cartn_y'},
                    $ref_atom_site->{$_}{'Cartn_z'} ] }
                ( $mid_atom_id, $up_atom_id, $side_atom_id );

        push @conformation_matrices,
             @{ bond_altering( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               $dihedral_angle_name,
                               $bond_angle_name,
                               $bond_name ) };
    }

    return \@conformation_matrices;
}

1;
