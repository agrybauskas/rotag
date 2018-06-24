package AlterMolecule;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( angle_bending
                     bond_positioning
                     bond_stretching
                     bond_torsion );

use Math::Trig;

use LinearAlgebra qw( mult_matrix_product
                      switch_ref_frame );
use Measure qw( bond_angle
                dihedral_angle );

# ----------------------- Molecule alteration matrices ------------------------ #

#
# Makes a rotational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_name - name of the dihedral angle.
# Output:
#     $rot_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined dihedral angle variables).
#

sub bond_torsion
{
    my ( $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $angle_name ) = @_;

    # Rotation matrix around the bond.
    my @rot_matrix =
        ( [ "cos(\$${angle_name})", "-sin(\$${angle_name})", 0, 0 ],
          [ "sin(\$${angle_name})",  "cos(\$${angle_name})", 0, 0 ],
          [ 0, 0, 1, 0 ],
          [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
        mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'global' ) },
                               \@rot_matrix,
                               @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'local' ) } ] );

    return $rot_matrix;
}

#
# Makes a translational transformation matrix for the bond of interest.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $length_name - name of the bond length variable.
# Output:
#     $transl_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined bond length variable).
#

sub bond_stretching
{
    my ( $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $length_name ) = @_;

    # Translation of the coordinates of the bond.
    my @transl_matrix = ( [ 1, 0, 0, 0 ],
                          [ 0, 1, 0, 0 ],
                          [ 0, 0, 1, "\$${length_name}" ],
                          [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $transl_matrix =
        mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'global' ) },
                               \@transl_matrix,
                               @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'local' ) } ] );

    return $transl_matrix;
}

#
# Makes a transformation matrix for changing angles between two bonds.
# Input:
#     ${mid,up,side}_atom_coord - Cartesian coordinates in array form that define
#     user-selected mid, up, side atoms.
#     $angle_name_x - name of the bond angle that will be rotated in x-axis.
#     $angle_name_y - name of the bond angle that will be rotated in y-axis.
# Output:
#     $rot_matrix - matrix defining coordinates in symbolic mathematical form
#     (with undefined angle variables).
#

sub angle_bending
{
    my ( $mid_atom_coord,
         $up_atom_coord,
         $side_atom_coord,
         $angle_name_x,
         $angle_name_y, ) = @_;

    # Bond angle matrices that rotates along x and y axes.
    my @rot_matrix_x =
        ( [ 1, 0, 0, 0 ],
          [ 0, "cos(\$${angle_name_x})", "-sin(\$${angle_name_x})", 0 ],
          [ 0, "sin(\$${angle_name_x})", "cos(\$${angle_name_x})", 0 ],
          [ 0, 0, 0, 1 ] );
    my @rot_matrix_y =
        ( [ "cos(\$${angle_name_y})", 0, "sin(\$${angle_name_y})", 0 ],
          [ 0, 1, 0, 0 ],
          [ "-sin(\$${angle_name_y})", 0, "cos(\$${angle_name_y})", 0 ],
          [ 0, 0, 0, 1 ] );

    # Multiplying multiple matrices to get a final form.
    my $rot_matrix =
        mult_matrix_product( [ @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'global' ) },
                               \@rot_matrix_y,
                               \@rot_matrix_x,
                               @{ switch_ref_frame( $mid_atom_coord,
                                                    $up_atom_coord,
                                                    $side_atom_coord,
                                                    'local' ) } ] );

    return $rot_matrix;
}


#            sp3                       sp2               sp
#
#            Up(2)                     Up(2)             Up(2)
# z          |                         |                 |
# |_y      Middle(1) __ Right(3)     Middle(1)         Middle(1)
# /         / \                       / \                |
# x    Left(4) Back(5)           Left(4) Right(3)        Down(3)
#

sub bond_positioning
{
    my ( %args ) = @_;
    my ( $mid_atom_coord,
         $up_atom_coord,
         $right_atom_coord,
         $left_atom_coord,
         $back_atom_coord,
         $add,
         $bond_length ) = (
        $args{'mid_atom_coord'},
        $args{'up_atom_coord'},
        $args{'right_atom_coord'},
        $args{'left_atom_coord'},
        $args{'back_atom_coord'},
        $args{'add'},
        $args{'bond_length'}
    );

    # Theoretically optimal angles should be so, the distance of
    # atoms would be furthest from each other. This strategy could
    # be achieved by imagining each bond as vector of the length 1.
    # Then the sum of all vectors would be 0, if angles are equal.
    #
    #                        ->  ->  ->  ->
    #                        A + B + C + D = 0
    #
    # In that case, the square of sum also should be equal 0.
    #
    #                        ->  ->  ->  ->
    #                      ( A + B + C + D ) ^ 2 = 0
    #
    #       ->  ->      ->  ->      ->  ->      ->  ->      ->  ->
    #   2 * A * B + 2 * A * C + 2 * A * D + 2 * B * C + 2 * B * D +
    #       ->  ->   ->      ->      ->      ->
    # + 2 * C + D  + A ^ 2 + B ^ 2 + C ^ 2 + D ^ 2 = 0
    #
    # Because length of each vector is equal to 1, then dot product
    # is equal to cos(alpha), where all angles between bonds are
    # equal. If there were no given angles, alpha should be 109.5
    # degrees.
    #
    #                   alpha = arccos( - 1 / 3 )
    #
    # However, there is a restriction of given angle. And the
    # calculation changes:
    #
    #        alpha = arccos( ( - 4 - 2 * cos( beta )
    #                              - 2 * cos( gamma )
    #                              - 2 * cos( delta ) ) / 6 )
    #
    # where beta is the given angle.

    # Calculates all bond angles between atoms connected to the
    # middle atom.
    my $up_mid_right_angle =
        bond_angle( [ $up_atom_coord,
                      $mid_atom_coord,
                      $right_atom_coord ] );
    my $up_mid_left_angle =
        bond_angle( [ $up_atom_coord,
                      $mid_atom_coord,
                      $left_atom_coord ] );
    my $right_mid_left_angle =
        bond_angle( [ $up_atom_coord,
                      $mid_atom_coord,
                      $left_atom_coord ] );

    # Calculates what common angle between atom and rest of the
    # atoms should be.
    my $bond_angle =
        acos( ( - 4
                - 2 * cos( $up_mid_right_angle )
                - 2 * cos( $up_mid_left_angle )
                - 2 * cos( $right_mid_left_angle ) )
              / 6 );

    # Determines dihedral angle between left and right atoms. Then
    # splits rest of the 2 * pi angle into two equal parts.
    my $dihedral_angle =
        dihedral_angle( [ $left_atom_coord,
                          $up_atom_coord,
                          $mid_atom_coord,
                          $right_atom_coord ] );
    if( abs( $dihedral_angle ) < ( 3 * pi() / 4 ) ) {
        if( $dihedral_angle < 0 ) {
            $dihedral_angle = ( 2 * pi() + $dihedral_angle ) / 2;
        } else {
            $dihedral_angle = - ( 2 * pi() - $dihedral_angle ) / 2;
        }
    } else {
        if( $dihedral_angle < 0 ) {
            $dihedral_angle = $dihedral_angle / 2;
        } else {
            $dihedral_angle = - $dihedral_angle / 2;
        }
    }

    # Places atom/moiety according to previously calculated angles.
    my ( $transf_matrix ) =
        @{ switch_ref_frame(
               $mid_atom_coord,
               $up_atom_coord,
               $left_atom_coord,
               'global' ) };

    $transf_matrix =
        mult_matrix_product(
            [ $transf_matrix,
              [ [ $bond_length
                  * cos( pi() / 2 - $dihedral_angle )
                  * sin( $bond_angle ) ],
                [ $bond_length
                  * sin( pi() / 2 - $dihedral_angle )
                  * sin( $bond_angle ) ],
                [ $bond_length
                  * cos( $bond_angle ) ],
                [ 1 ] ] ] );

    return $transf_matrix;
}

1;
