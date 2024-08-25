package Grid;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_box
                     grid_box
                     identify_neighbour_cells );

use List::Util qw( max
                   min );
use Carp qw( confess );
use ForceField::Parameters;
use PDBxParser qw( extract );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------- Grid-related functions --------------------------- #

#
# Given the cartesian coordinates (x, y, z) of atoms, function returns the
# dimensions of smallest possible box that contains all atoms.
# Input:
#     @atom_coord - list of atom coordinates in x, y, z array form.
# Output:
#     coordinates of min and max x, y, z box boundaries in which all given atoms
#     are contained.
#

sub create_box
{
    my ( $atom_coord ) = @_;

    my @atom_coord_x = map { $_->[0] } @{ $atom_coord };
    my @atom_coord_y = map { $_->[1] } @{ $atom_coord };
    my @atom_coord_z = map { $_->[2] } @{ $atom_coord };

    if( ! ( defined $atom_coord_x[0] &&
            defined $atom_coord_y[0] &&
            defined $atom_coord_z[0] ) ) {
        confess 'not all x, y and z coordinates are defined';
    }

    # Directions are adapted to right-handed Cartesian coordinate system.
    # Looking for leftmost and rightmost coordinates of X-axis.
    my $min_coord_x  = min( @atom_coord_x );
    my $max_coord_x  = max( @atom_coord_x );

    # Looking for most backward and forward coordinates of Y-axis.
    my $min_coord_y  = min( @atom_coord_y );
    my $max_coord_y  = max( @atom_coord_y );

    # Looking for downmost and upmost coordinates of Z-axis.
    my $min_coord_z  = min( @atom_coord_z );
    my $max_coord_z  = max( @atom_coord_z );

    # Coordinates of minimum bounding box that contains all given atoms.
    return [ $min_coord_x, $max_coord_x,
	     $min_coord_y, $max_coord_y,
	     $min_coord_z, $max_coord_z ];
}

#
# Divides atoms into grid box of given edge length.
# Input:
#     $atom_site - special atom data structure (see PDBxParser.pm);
#     $edge_length - edge length of the cell inside grid box.
# Output:
#     %grid_box - hash where key is string representing cell id and value -
#     atom id.
#

sub grid_box
{
    my ( $parameters, $atom_site, $edge_length, $atom_ids, $options ) = @_;
    my ( $attributes ) = ( $options->{'attributes'} );

    my $edge_length_connection =
        $parameters->{'_[local]_constants'}{'edge_length_connection'};

    $attributes //= [ 'id' ];

    # Default value for edge length is two times greater than the largest
    # covalent radius.
    $edge_length //= $edge_length_connection;

    # Determines boundary box around all atoms.
    my $atom_data =
        extract( $atom_site,
                 { 'data' => [ 'id', 'Cartn_x', 'Cartn_y', 'Cartn_z' ] } );
    my @atom_coordinates = map { [ $_->[1], $_->[2], $_->[3] ] } @{ $atom_data };
    my $boundary_box = create_box( \@atom_coordinates );

    # Creates box with cells with edge length of given variable in angstroms.
    my %grid_box;
    my %atom_cell_pos;
    my $cell_index_x;
    my $cell_index_y;
    my $cell_index_z;

    # Iterates through atoms and determines in which cell these atoms are.
    foreach my $atom_coord ( @{ $atom_data } ) {
        $cell_index_x =
            int( ( $atom_coord->[1] - $boundary_box->[0] ) / $edge_length ) + 1;
        $cell_index_y =
            int( ( $atom_coord->[2] - $boundary_box->[2] ) / $edge_length ) + 1;
        $cell_index_z =
            int( ( $atom_coord->[3] - $boundary_box->[4] ) / $edge_length ) + 1;

        # Checks if hash keys already  exist.
        if( exists $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} ) {
            push @{ $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} },
                $atom_coord->[0];
        } else {
            $grid_box{"$cell_index_x,$cell_index_y,$cell_index_z"} =
                [ $atom_coord->[0] ];
        }

        # Identifies what cells do atoms occupy, if they are in the list.
        if( ! $atom_ids ) { next; }

        for my $atom_id ( @{ $atom_ids } ) {
        if( $atom_coord->[0] eq $atom_id ) {
            push @{ $atom_cell_pos{"$cell_index_x,$cell_index_y,$cell_index_z"}},
                 join q(,), map { $atom_site->{$atom_coord->[0]}{$_} }
                 @{ $attributes };
        } }
    }

    return \%grid_box, \%atom_cell_pos;
}

#
# Identifies neighbouring cells of selected cell.
# Input:
#     $grid_box - grid box generated by function grid_box;
#     $specified_cells - list of cell ids.
# Output:
#     %neighbour_cells - hash of neighbouring cells with atom ids inside
#     (data structure is the same as $grid_box).
#

sub identify_neighbour_cells
{
    my ( $grid_box, $specified_cells ) = @_;

    $specified_cells //= $grid_box; # Checks all cells.

    # Checks for neighbouring cells for each cell.
    my %neighbour_cells;
    foreach my $cell ( keys %{ $specified_cells } ) {
        my @cell_idxs = split /,/sxm, $cell;

        # $i represents x, $j - y, $k - z coordinates.
        for my $i ( ( $cell_idxs[0] - 1..$cell_idxs[0] + 1 ) ) {
        for my $j ( ( $cell_idxs[1] - 1..$cell_idxs[1] + 1 ) ) {
        for my $k ( ( $cell_idxs[2] - 1..$cell_idxs[2] + 1 ) ) {
            if( exists $grid_box->{"$i,$j,$k"} ) {
                push @{ $neighbour_cells{"$cell"} },
                     @{ $grid_box->{"$i,$j,$k"} };
        } } } }
    }

    return \%neighbour_cells;
}

1;