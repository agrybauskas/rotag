package forseeAtomPos;

use strict;
use warnings;

# ------------------------------- Linear algebra ------------------------------ #

#
# Description. Creates local reference frame for any three given atoms.
# Input: 3x3 array, where each row represents x, y, z coordinates of one of three
#        atoms.
# Output: 3x3 array, where each row represents x, y, z coordinates of 
#         corresponding x-axis, y-axis, z-axis of local frame of reference.
#

sub createRefFrame{
    my ($up_atom_x,   $up_atom_y,   $up_atom_z,
        $side_atom_x, $side_atom_y, $side_atom_z,
        $mid_atom_x,  $mid_atom_y,  $mid_atom_z) = @_;
    
    my @local_ref;
    
    # Let local z-axis be colinear to bond between mid and up atoms.
    $local_ref_frame[2][0] = $side_atom_x - $up_atom_x;
    $local_ref_frame[2][1] = $side_atom_y - $up_atom_y;
    $local_ref_frame[2][2] = $side_atom_z - $up_atom_z;
    
    # Let local x-axis be perpendicular to bonds between mid, up and mid, side
    # atoms.
    $local_ref_frame[0][0] = ($mid_atom_y - $up_atom_y) * $local_ref_frame[2][2]
        - ($mid_atom_z - $up_atom_z) * $local_ref_frame[2][1];
    $local_ref_frame[0][1] = - ($mid_atom_x - $up_atom_x) * $local_ref_frame[2][2]
        + ($mid_atom_z - $up_atom_z) * $local_ref_frame[2][0];
    $local_ref_frame[0][2] = ($mid_atom_x - $up_atom_x) * $local_ref_frame[2][1]
        - ($mid_atom_y - $up_atom_y) * $local_ref_frame[2][0];
    
    # Let local y-axis be in the same plane as mid-up and mid-side bonds.
    $local_ref_frame[1][0] = $local_ref_frame[2][1] * $local_ref_frame[0][2]
        - $local_ref_frame[2][2] * $local_ref_frame[0][1];
    $local_ref_frame[1][1] =-$local_ref_frame[2][0]*$local_ref_frame[0][2]
        + $local_ref_frame[2][2] * $local_ref_frame[0][0];
    $local_ref_frame[1][2] = $local_ref_frame[2][0] * $local_ref_frame[0][1]
        - $local_ref_frame[2][1] * $local_ref_frame[0][0];

    return @local_ref_frame;
}

1;
