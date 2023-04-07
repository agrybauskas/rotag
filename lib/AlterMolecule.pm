package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending
                     bond_altering
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
#     $angle_name - name of the dihedral angle.
# Output:
#     $rot_matrix - matrix defining coordinates in analytical form.
#

sub bond_torsion
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $angle_name ) = @_;

    # Rotation matrix around the bond.
    my $rot_matrix =
        Symbolic->new(
            { 'symbols' => [ $angle_name ],
              'matrix' => sub { my ( $svar ) = @_;
                                return [ [ cos( $svar ),-sin( $svar ), 0, 0 ],
                                         [ sin( $svar ), cos( $svar ), 0, 0 ],
                                         [ 0, 0, 1, 0 ],
                                         [ 0, 0, 0, 1 ], ]; } } );

    my @rot_matrix =
        ( @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'global' ) },
          $rot_matrix,
          @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'local' ) }, );

    return \@rot_matrix;
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

#
# Makes a transformation matrix for changing bond lengths, angles and dihedral
# angles.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms;
#     $dihedral_angle_name - name of the dihedral angle that will be rotated in
#     z-axis;
#     $bond_angle_name - name of the bond angle that will be rotated in y-axis;
#     $bond_name - name of the bond that will be translated along z-axis.
# Output:
#     $altering_matrix - matrix defining coordinates in analytical form.
#

sub bond_altering
{
    my ( $parameters,
         $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $dihedral_angle_name,
         $bond_angle_name,
         $bond_name ) = @_;

    # Rotation matrix around the bond.
    # TODO: a_{i-1}, - d_i * sin(alpha_{i-1}, d_i * cos(alpha_{i-1}).
    my $bond_altering_matrix =
        Symbolic->new(
            { 'symbols' => [ $dihedral_angle_name, $bond_angle_name, $bond_name ],
              'matrix' =>
                  sub { my ( $svar1, $svar2, $svar3 ) = @_;
                        return [ [ cos( $svar1 ), -sin( $svar1 ), 0, $up_atom_coord->[2] - $mid_atom_coord->[2] ],
                                 [ sin( $svar1 ) * cos( $svar2 ), cos( $svar1 ) * cos( $svar2 ), -sin( $svar2 ), 0 ],
                                 [ sin( $svar1 ) * sin( $svar2 ), cos( $svar1 ) * sin( $svar2 ),  cos( $svar2 ), 0 ],
                                 [ 0, 0, 0, 1 ], ]; } } );

    my @bond_altering_matrix =
        ( @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'global' ) },
          $bond_altering_matrix,
          @{ switch_ref_frame( $parameters,
                               $mid_atom_coord,
                               $up_atom_coord,
                               $side_atom_coord,
                               'local' ) }, );

    return \@bond_altering_matrix;
}

1;
