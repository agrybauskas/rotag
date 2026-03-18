package DBScan;

use strict;
use warnings;

use List::Util qw( max
                   sum );

use DBScan;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $points ) = @_;
    my $self = {};
    my $id = 1;
    for my $point ( @{ $points } ) {
        $self->{$id}{'id'} = $id;
        $self->{$id}{'label'} = undef;
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
# Adds points.
# Input:
#     $coord - coordinates;
#     $id - point id;
# Output:
#     appends points to the data structure.
#

sub add_point
{
    my ( $self, $coord, $id ) = @_;
    $id //= max @{ $self->get_ids() } + 1;
    $self->{$id}{'id'} = $id;
    $self->{$id}{'label'} = undef;
    $self->{$id}{'coord'} = $coord;
    return;
}

#
# Performs DBSCAN clusterization for given data points.
# Input:
#     $radius - distance between points to be included to the cluster;
#     $options->{'min_points'} - minimum points required for them to be in the
#     cluster;
# Output:
#     labels points by cluster name.
#

sub dbscan
{
    my ( $self, $radius, $options ) = @_;
    my ( $min_points ) = ( $options->{'min_points'} );
    $min_points //= 1;
    my $label = 0;
    for my $id ( @{ $self->get_ids() } ) {
        next if defined $self->{$id}{'label'} && $self->{$id}{'label'} > 0;

        my $neighbours = range_query( $self, $self->{$id}, $radius );

        if( scalar @{ $neighbours } < $min_points ) {
            $self->{$id}{'label'} = -1;
        }

        $label++;
        $self->{$id}{'label'} = $label;

        for my $neighbour_id ( @{ $neighbours }) {
            next if ! defined $self->{$neighbour_id}{'label'};

            if( $self->{$neighbour_id}{'label'} < 0 ) {
                $self->{$neighbour_id}{'label'} = $label;
            }

            if( $self->{$neighbour_id}{'label'} > 0 ) {
                $self->{$id}{'label'} = $self->{$neighbour_id}{'label'};
                $label--;  # Resets the label as it was not used.
            }
        }
    }

    return;
}

# -------------------------------- Functions ---------------------------------- #

#
# Detects neighbouring points by defined search radius.
# Input:
#     $points - data structure created by DBScan library containing
#     coordinate values, cluster names and ids;
#     $point - a point from $points data structure which neighbours are searched
#     for;
# Output:
#     @neighbour_ids -- return neighbouring point ids.
#

sub range_query
{
    my ( $points, $point, $radius ) = @_;
    my @neighbour_ids = ();
    my $radius_squared = $radius ** 2;
    for my $point_id ( @{ $points->get_ids() } ) {
        next if $point_id eq $point->{'id'};

        my $dimension = scalar @{ $point->{'coord'} };

        my $distance_squared =
            sum( map { ( $points->{$point_id}{'coord'}[$_] -
                         $point->{'coord'}[$_] ) ** 2 }
                     ( 0..$dimension-1 ) );

        next if $distance_squared > $radius_squared;

        push @neighbour_ids, $point_id;
    }
    return \@neighbour_ids;
}

1;
