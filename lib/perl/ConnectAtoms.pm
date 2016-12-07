package ConnectAtoms;

use strict;
use warnings;

use List::Util qw(min max);
use POSIX qw(ceil);
use List::MoreUtils qw(zip);
use Data::Dumper;

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
# Input   (1 arg): atom coordinates in x, y, z form.
# Output (6 args): coordinates of min and max x, y, z box boundaries in which
#                  all given atoms are contained.
#

sub create_box
{
    my @atom_coord = @_;

    my @atom_coord_x = map { $_->[0] } @atom_coord;
    my @atom_coord_y = map { $_->[1] } @atom_coord;
    my @atom_coord_z = map { $_->[2] } @atom_coord;

    # Directions are adapted to right-handed Cartesian coordinate system.
    # Looking for leftmost and rightmost coordinates of X-axis.
    my $most_left_x_coord     = min( @atom_coord_x );
    my $most_right_x_coord    = max( @atom_coord_x );

    # Looking for most backward and forward coordinates of Y-axis.
    my $most_backward_y_coord = min( @atom_coord_y );
    my $most_forward_y_coord  = max( @atom_coord_y );

    # Looking for downmost and upmost coordinates of Z-axis.
    my $most_down_z_coord     = min( @atom_coord_z );
    my $most_up_z_coord       = max( @atom_coord_z );

    # Coordinates of minimum bounding box that contains all given atoms.
    return $most_left_x_coord,     $most_right_x_coord,
           $most_backward_y_coord, $most_forward_y_coord,
           $most_down_z_coord,     $most_up_z_coord;
}

#
# Divides box into grid of cubes that has length of the desired bond. If box
# is not perfectly divisible, then the boundaries are extended accordingly.
# Input (2 args): bond length in angstroms, coordinates of atoms (x, y, z).
#

sub connect_atoms
{
    my $bond_length = shift;
    my @all_atom_coord  = @_;

    # Create smallest box that contain all atoms.
    my @boundary_box = create_box( @all_atom_coord );

    # Divide box into cells, edges of which has length of bond.
    my $cell_number_x =
	ceil( ( $boundary_box[1] - $boundary_box[0] ) / $bond_length );
    my $cell_number_y =
	ceil( ( $boundary_box[3] - $boundary_box[2] ) / $bond_length );
    my $cell_number_z =
	ceil( ( $boundary_box[5] - $boundary_box[4] ) / $bond_length );

    my @grid_box;

    foreach my $idx_x ( ( 1..$cell_number_x ) ) {
	foreach my $idx_y ( ( 1..$cell_number_y ) ) {
	    foreach my $idx_z ( ( 1..$cell_number_z ) ) {
		push( @grid_box, [ $idx_x, $idx_y, $idx_z ] => [] );
	    }
	}
    }

    # Assign atoms to cells in grid_box.
    foreach my $atom_coord ( @all_atom_coord ) {
	foreach my $cell ( @grid_box ) {
	}
    }

}

1;
