package ConnectAtoms;

use List::Util qw(min max);

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
#

sub create_box
{
    my @atom_coordinates = @_;

    # Directions are adapted to right-handed Cartesian coordinate system.
    # Looking for leftmost and rightmost coordinates of X-axis.
    my $most_left_x_coord     = min( $atom_coordinates[0] );
    my $most_right_x_coord    = max( $atom_coordinates[0] );

    # Looking for most backward and forward coordinates of Y-axis.
    my $most_backward_y_coord = min( $atom_coordinates[1] );
    my $most_forward_y_coord  = max( $atom_coordinates[1] );

    # Looking for downmost and upmost coordinates of Z-axis.
    my $most_down_z_coord     = min( $atom_coordinates[2] );
    my $most_up_z_coord       = max( $atom_coordinates[2] );

    # Coordinates of minimum bounding box that contains all given atoms.
    my @box_coordinates;

    $box_coordinates[0] =
	[ $most_left_x_coord,  $most_backward_y_coord, $most_down_z_coord ];
    $box_coordinates[1] =
	[ $most_left_x_coord,  $most_backward_y_coord, $most_up_z_coord   ];
    $box_coordinates[2] =
	[ $most_left_x_coord,  $most_forward_y_coord,  $most_up_z_coord   ];
    $box_coordinates[3] =
    	[ $most_left_x_coord,  $most_forward_y_coord,  $most_down_z_coord ];
    $box_coordinates[4] =
    	[ $most_right_x_coord, $most_backward_y_coord, $most_down_z_coord ];
    $box_coordinates[5] =
    	[ $most_right_x_coord, $most_backward_y_coord, $most_up_z_coord   ];
    $box_coordinates[6] =
    	[ $most_right_x_coord, $most_forward_y_coord,  $most_up_z_coord   ];
    $box_coordinates[7] =
    	[ $most_right_x_coord, $most_forward_y_coord,  $most_down_z_coord ];

    print Dumper \@box_coordinates;

    return \@box_coordinates;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
#

sub divide_box
{

}

1;
