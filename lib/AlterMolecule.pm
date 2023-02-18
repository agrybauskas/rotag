package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending
                     bond_stretching
                     bond_torsion );

use LinearAlgebra qw( switch_ref_frame );
use Symbolic;
use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------- Molecule alteration matrices ------------------------ #

#
# Makes a rotational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms;
#     $angle_name_z - name of the dihedral angle.
# Output:
#     $rot_matrix_z - matrix defining coordinates in analytical form.
#

sub bond_torsion
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $angle_name_z,
         $angle_name_x ) = @_;

    # Rotation matrix around the bond.
    my $rot_matrix_z =
        Symbolic->new(
            { 'symbols' => [ $angle_name_z ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ cos( $svar ),-sin( $svar ), 0, 0 ],
                                         [ sin( $svar ), cos( $svar ), 0, 0 ],
                                         [ 0, 0, 1, 0 ],
                                         [ 0, 0, 0, 1 ], ]; } } );
    # x-axis rotation is used in combination with bond angle transformation,
    # because it changes the atom positions that are used to calculate the
    # frame of reference.
    my $rot_matrix_x =
        Symbolic->new(
            { 'symbols' => [ $angle_name_x ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ 1, 0, 0, 0 ],
                                         [ 0, cos( $svar ), -sin( $svar ), 0 ],
                                         [ 0, sin( $svar ), cos( $svar ), 0 ],
                                         [ 0, 0, 0, 1 ], ]; } } );

    my @rot_matrix_z =
        ( @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'global' ) },
          $rot_matrix_z,
          ( defined $angle_name_x ? $rot_matrix_x : ()),
          @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'local' ) }, );

    return \@rot_matrix_z;
}

#
# Makes a translational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms;
#     $length_name - name of the bond length variable.
# Output:
#     $transl_matrix - matrix defining coordinates in analytical form.
#

sub bond_stretching
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $length_name ) = @_;

    # Translation of the coordinates of the bond.
    my $transl_matrix =
        Symbolic->new(
            { 'symbols' => [ $length_name ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ 1, 0, 0, 0 ],
                                         [ 0, 1, 0, 0 ],
                                         [ 0, 0, 1, $svar ],
                                         [ 0, 0, 0, 1 ], ]; } } );

    # Multiplying multiple matrices to get a final form.
    my @transl_matrix =
        ( @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'global' ) },
          $transl_matrix,
          @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'local' ) }, );

    return \@transl_matrix;
}

#
# Makes a transformation matrix for changing angles between two bonds.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms;
#     $angle_name_x - name of the bond angle that will be rotated in x-axis;
#     $angle_name_y - name of the bond angle that will be rotated in y-axis.
# Output:
#     $rot_matrix - matrix defining coordinates in analytical form.
#

sub angle_bending
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $angle_name_x,
         $angle_name_y ) = @_;

    # Bond angle matrices that rotates along x and y axes.
    my $rot_matrix_x =
        Symbolic->new(
            { 'symbols' => [ $angle_name_x ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ 1, 0, 0, 0 ],
                                         [ 0, cos( $svar ), -sin( $svar ), 0 ],
                                         [ 0, sin( $svar ), cos( $svar ), 0 ],
                                         [ 0, 0, 0, 1 ], ]; } } );
    my $rot_matrix_y =
        Symbolic->new(
            { 'symbols' => [ $angle_name_y ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ cos( $svar ), 0, sin( $svar ), 0 ],
                                         [ 0, 1, 0, 0 ],
                                         [ -sin( $svar ), 0, cos( $svar ), 0 ],
                                         [ 0, 0, 0, 1 ], ]; } } );

    # Multiplying multiple matrices to get a final form.
    my @rot_matrix =
        ( @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'global' ) },
          $rot_matrix_y,
          $rot_matrix_x,
          @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'local' ) }, );

    return \@rot_matrix;
}

1;
