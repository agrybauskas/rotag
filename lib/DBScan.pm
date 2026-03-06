package DBScan;

use strict;
use warnings;

use DBScan;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $points ) = @_;
    my $self = {};
    my $id = 1;
    for my $point ( @{ $points } ) {
        $self->{$id}{'id'} = $id;
        $self->{$id}{'label'} = -1;
        $self->{$id}{'coord'} = $point;
        $id++;
    }
    return bless $self, $class;
}

# --------------------------------- Methods ----------------------------------- #

#
# Performs DBSCAN clusterization for given data points.
# Input:
#     $points - data structure that stores both point coordination and label;
#     $radius - distance between points to be included to the cluster;
#     $min_points - minimum points required by the cluster;
#     $distance_func - dinstance function to be calculated by;
# Output:
#     labels points by cluster name.
#

sub dbscan
{
    my ( $self, $radius, $options ) = @_;
    my ( $min_points, $distance_func ) =
        ( $options->{'min_points'}, $options->{'distance_func'} );
    $min_points //= 1;
    for my $id ( sort keys %{ $self } ) {
        my $point = $self->{$id};

        next if $point->{'label'} > 0;

        my $neighbours = range_query( $self, $point, $distance_func, $radius );
    }
}

# -------------------------------- Functions ---------------------------------- #

sub range_query
{
    my ( $points, $point, $distance_func, $radius ) = @_;
}

1;
