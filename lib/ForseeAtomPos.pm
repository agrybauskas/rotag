package ForseeAtomPos;

use strict;
use warnings;

# ------------------------------ PDBx/mmCIF parser ---------------------------- #

sub cif_to_hash{
    my %cif_hash;
    my $read_atom_site = 0;  
    my @cif_table;  
    my $cif_line_counter = 0;

    foreach( @_ ){
        if( $_ =~ /_atom_site.(.+)\n$/ ){
            $cif_hash{ $_ }{ $cif_line_counter } = "";
            $read_atom_site = 1;
            $cif_line_counter += 1;
        } elsif( $read_atom_site == 1 ){
            push( @cif_table, split( /\s+/, $_ ) );
        } elsif( $read_atom_site == 1 && $_ =~ /[^_]/ ){
            last;
        }
    }

    
}

# ------------------------------- Linear algebra ------------------------------ #

#
# Constants
#

my $PI = 4 * atan2( 1, 1 );
my $EPSILON = 1.0 / ( 2 ** 52 ); # Machine accuracy for 64-bit floating point
                                 # numbers.

#
# Description. Creates local reference frame for any three given atoms.
# Input: 3x3 array, where each row represents x, y, z coordinates of one of three
#        atoms.
# Output: 3x3 array, where each row represents x, y, z coordinates of 
#         corresponding x-axis, y-axis, z-axis of local frame of reference.
#

sub create_ref_frame{
    my ( $mid_atom_x,  $mid_atom_y,  $mid_atom_z,
         $up_atom_x,   $up_atom_y,   $up_atom_z,
         $side_atom_x, $side_atom_y, $side_atom_z ) = @_;

    my @local_ref_frame;

    # Let local z-axis be colinear to bond between mid and up atoms.
    $local_ref_frame[2][0] = $up_atom_x - $mid_atom_x;
    $local_ref_frame[2][1] = $up_atom_y - $mid_atom_y;
    $local_ref_frame[2][2] = $up_atom_z - $mid_atom_z;

    # Let local x-axis be perpendicular to bonds between mid, up and mid, side
    # atoms.
    $local_ref_frame[0][0] = 
        ( $side_atom_y - $mid_atom_y ) * $local_ref_frame[2][2]
      - ( $side_atom_z - $mid_atom_z ) * $local_ref_frame[2][1];
    $local_ref_frame[0][1] = 
      - ( $side_atom_x - $mid_atom_x ) * $local_ref_frame[2][2]
      + ( $side_atom_z - $mid_atom_z ) * $local_ref_frame[2][0];
    $local_ref_frame[0][2] = 
        ( $side_atom_x - $mid_atom_x ) * $local_ref_frame[2][1]
      - ( $side_atom_y - $mid_atom_y ) * $local_ref_frame[2][0];

    # Let local y-axis be in the same plane as mid-up and mid-side bonds.
    $local_ref_frame[1][0] = 
        $local_ref_frame[2][1] * $local_ref_frame[0][2]
      - $local_ref_frame[2][2] * $local_ref_frame[0][1];
    $local_ref_frame[1][1] =
      - $local_ref_frame[2][0] * $local_ref_frame[0][2]
      + $local_ref_frame[2][2] * $local_ref_frame[0][0];
    $local_ref_frame[1][2] = 
        $local_ref_frame[2][0] * $local_ref_frame[0][1]
      - $local_ref_frame[2][1] * $local_ref_frame[0][0];

    return @local_ref_frame;
}

#
# Description. Function calculates Euler rotational angles (alpha, beta, gamma)
# that are used to transform global reference frame to chosen one.
# Input: 3x3 array, where each row represents vectors of x-axis, y-axis and 
#        z-axis accordingly.
# Output: alpha, beta, gamma angles in radians.
#

sub find_euler_angles{
    my ( $mid_atom_x,  $mid_atom_y,  $mid_atom_z,
         $up_atom_x,   $up_atom_y,   $up_atom_z,
         $side_atom_x, $side_atom_y, $side_atom_z ) = @_;

    my $alpha_rad;
    my $beta_rad;
    my $gamma_rad;

    my $z_axis_in_xy_plane;

    my @local_ref_frame = 
        create_ref_frame( $mid_atom_x,  $mid_atom_y,   $mid_atom_z,
                          $up_atom_x,   $up_atom_y,    $up_atom_z,
                          $side_atom_x, $side_atom_y,  $side_atom_z );

    # Projects local z-axis to global xy-plane.
    $z_axis_in_xy_plane = 
        sqrt( $local_ref_frame[2][0] * $local_ref_frame[2][0] 
            + $local_ref_frame[2][1] * $local_ref_frame[2][1] );

    if( $z_axis_in_xy_plane > $EPSILON ){
        $alpha_rad = 
            atan2( $local_ref_frame[1][0] * $local_ref_frame[2][1]
                 - $local_ref_frame[1][1] * $local_ref_frame[2][0],
                   $local_ref_frame[0][0] * $local_ref_frame[2][1]
                 - $local_ref_frame[0][1] * $local_ref_frame[2][0] );
        $beta_rad = atan2( $z_axis_in_xy_plane, $local_ref_frame[2][2] );
        $gamma_rad = - atan2( - $local_ref_frame[2][0], $local_ref_frame[2][1] );
    } else {
        $alpha_rad = 0.;
        $beta_rad = ( $local_ref_frame[2][2] > 0. ) ? 0. : $PI;
        $gamma_rad = - atan2( $local_ref_frame[0][1], $local_ref_frame[0][0] );
    }

    return $alpha_rad, $beta_rad, $gamma_rad;
}

# --------------------- Computer algebra software wrappers -------------------- #

#
# Because Perl (v5.14.2) is not capable of performing symbolic algebra these
# functions act as bridges/wrappers between Perl and programs or modules that
# can perform symbolic computations, such as Maxima, GNU Octave, GiNaC 
# (C++ package) and etc.
#
# Example of rotation along z-axis by chi angle in radians:
#
#      / cos(chi) -sin(chi) 0 \   / x \   / x * cos(chi) + y * sin(chi) \
#      | sin(chi)  cos(chi) 0 | * | y | = | x * sin(chi) + y * cos(chi) |
#      \    0         0     1 /   \ z /   \              0              /
# 
#
# Description. A wrapper function for Maxima 5.24.0. Takes argument from amino 
# acid model function and performs symbolic matrix multiplications with unknown 
# variables, such as bond dihedral chi angles, and simplifies the expression.
# Input: List of arrays representing matrices and strings indicating symbolic 
#        variables. The input comes from output of amino acid model functions. 
# Output: Simplified version of list of arrays representing matrices and 
#         strings indicating symbolic variables.
#

sub simplify_maxima{}

1;
