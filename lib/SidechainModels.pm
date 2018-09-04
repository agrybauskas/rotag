package SidechainModels;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotation_only );

use List::MoreUtils qw( uniq );

use AlterMolecule qw( bond_torsion );
use AtomProperties qw( sort_atom_names );
use LinearAlgebra qw( mult_matrix_product
                      reshape );
use BondProperties qw( rotatable_bonds );
use PDBxParser qw( filter
                   filter_by_unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------- Idealistic sidechain models ------------------------- #

#
# Model that uses only rotation around single bonds.
# Input:
#     $atom_site - atom data structure.
# Output:
#     adds conformation variable for each atom as list of transformation
#     matrices.
#

sub rotation_only
{
    my ( $atom_site ) = @_;

    my %atom_site = %{ $atom_site }; # Copy of $atom_site.

    # Determines all residue ids present in atom site.
    my @residue_unique_keys =
        uniq map { join q{,}, @{$_} }
              @{ filter( { 'atom_site' => $atom_site,
                           'data' => [ 'label_seq_id',
                                       'label_asym_id',
                                       'label_entity_id',
                                       'label_alt_id' ] } ) };

    # Iterates through target residues and their atom ids and assigns
    # conformational equations which can produce pseudo-atoms later.
    for my $residue_unique_key ( @residue_unique_keys ) {
        my $residue_site =
            filter_by_unique_residue_key( \%atom_site, $residue_unique_key );

        my $rotatable_bonds = rotatable_bonds( $residue_site );

        for my $atom_id ( keys %{ $residue_site }  ) {
            my @atom_coord = ( $atom_site{"$atom_id"}{'Cartn_x'},
                               $atom_site{"$atom_id"}{'Cartn_y'},
                               $atom_site{"$atom_id"}{'Cartn_z'}, );

            if( ! exists $rotatable_bonds->{$atom_id} ) { next; }

            my @transf_matrices; # Matrices for transforming atom coordinates.

            for my $angle_name ( sort { $a cmp $b }
                                 keys %{ $rotatable_bonds->{$atom_id} } ) {
                # First, checks if rotatable bond has fourth atom produce
                # dihedral angle. It is done by looking at atom connections: if
                # rotatable bond ends with terminal atom, then this bond is
                # excluded.
                my $up_atom_id = $rotatable_bonds->{$atom_id}{$angle_name}[1];
                if( scalar( @{ $residue_site->{$up_atom_id}
                                              {'connections'} } ) < 2 ){ next; }

                my $mid_atom_id = $rotatable_bonds->{$atom_id}{$angle_name}[0];
                my @mid_connections = # Excludes up atom.
                    grep { $_ ne $up_atom_id }
                    @{ $residue_site->{$mid_atom_id}{'connections'} };
                my @mid_connection_names = # Excludes up atom.
                    map { $residue_site->{$_}{'label_atom_id'} }
                    @mid_connections;
                my $side_atom_name =
                    sort_atom_names( \@mid_connection_names )->[0];
                my $side_atom_id =
                    filter( { 'atom_site' => $residue_site,
                              'include' =>
                            { 'label_atom_id' => [ $side_atom_name ] },
                              'data' => [ 'id' ],
                              'is_list' => 1 } )->[0];

                my $mid_atom_coord =
                    [ $residue_site->{$mid_atom_id}{'Cartn_x'},
                      $residue_site->{$mid_atom_id}{'Cartn_y'},
                      $residue_site->{$mid_atom_id}{'Cartn_z'} ];
                my $up_atom_coord =
                    [ $residue_site->{$up_atom_id}{'Cartn_x'},
                      $residue_site->{$up_atom_id}{'Cartn_y'},
                      $residue_site->{$up_atom_id}{'Cartn_z'} ];
                my $side_atom_coord =
                    [ $residue_site->{$side_atom_id}{'Cartn_x'},
                      $residue_site->{$side_atom_id}{'Cartn_y'},
                      $residue_site->{$side_atom_id}{'Cartn_z'} ];

            # Creates and appends matrices to a list of matrices that later
            # will be multiplied.
            push @transf_matrices,
                 @{ bond_torsion( $mid_atom_coord,
                                  $up_atom_coord,
                                  $side_atom_coord,
                                  $angle_name ) };
            }

            $atom_site->{$atom_id}{'conformation'} =
                mult_matrix_product(
                    [ @transf_matrices,
                      @{ reshape( [ @atom_coord, 1 ], [ 4, 1 ] ) } ] );
        }
    }

    return;
}

1;
