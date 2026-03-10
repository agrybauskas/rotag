package DBScan;

use strict;
use warnings;

use List::Util qw( sum );

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

# ----------------------------- Setters/Getters ------------------------------- #

sub get_ids
{
    my ( $self ) = @_;
    return [ sort { $a <=> $b } keys %{ $self } ];
}

sub get_point
{
    my ( $self, $id ) = @_;
    return $self->{$id};
}

# --------------------------------- Methods ----------------------------------- #

#
# Performs DBSCAN clusterization for given data points.
# Input:
#     $points - data structure that stores both point coordination and label;
#     $radius - distance between points to be included to the cluster;
#     $min_points - minimum points required by the cluster;
# Output:
#     labels points by cluster name.
#

sub dbscan
{
    my ( $self, $radius, $options ) = @_;
    my ( $min_points ) = ( $options->{'min_points'} );
    $min_points //= 1;
    my $label = 1;
    for my $id ( sort keys %{ $self } ) {
        my $point = $self->{$id};

        next if $point->{'label'} > 0;

        my $neighbours = range_query( $self, $point, $radius );
        for my $neighbour_id ( @{ $neighbours }) {
            if( $self->{$neighbour_id}{'label'} > 0 ) {
                $self->{$id}{'label'} = $self->{$neighbour_id}{'label'};
            } else {

            }
        }
    }
}

# -------------------------------- Functions ---------------------------------- #

sub range_query
{
    my ( $points, $point, $radius ) = @_;
    my @neighbour_ids = ();
    my $radius_squared = $radius ** 2;
    for my $point_id ( sort { $a <=> $b } keys %{ $points } ) {
        next if $point_id eq $point->{'id'};

        my $distance_squared =
            sum( map { $_ ** 2 } @{ $points->{$point_id}{'coord'} } );

        next if $distance_squared > $radius_squared;

        push @neighbour_ids, $point_id;
    }
    return \@neighbour_ids;
}

1;
