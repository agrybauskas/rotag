package DBScan;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw();

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------------- DBScan ---------------------------------- #

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
    my ( $points, $radius, $options ) = @_;
    my ( $min_points, $distance_func ) =
        ( $options->{'min_points'}, $options->{'distance_func'} );
    $min_points //= 1;
}

1;
