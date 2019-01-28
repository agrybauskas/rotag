package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles
                     sample_angles_qs_parsing );

use POSIX;

use Constants qw( $PI );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Sampling ---------------------------------- #

#
# Produces angle values that are separated by even intervals.
# Input:
#     $angle_ranges - boundary between which angles can be sampled.
#     $small_angle - smallest angle increment.
# Output:
#     @angles - sampled angles.
#
#

sub sample_angles
{
    my ( $angle_ranges, $small_angle ) = @_;

    my @angles;
    my $min_angle;
    my $max_angle;

    # Devides full circle (2*pi) into even intervals by $small_angle value.
    $small_angle = # Adjusts angle so, it could be devided evenly.
        2 * $PI / floor( 2 * $PI / $small_angle );
    my @small_angles =
        map { $_ * $small_angle } 0..( floor( 2 * $PI / $small_angle ) - 1 );

    # Iterates around the circle and adds evenly spaced angles, if they are
    # inside intervals ($angle_ranges).
    for my $angle ( @small_angles ) {
        # TODO: might speed up calculation by eliminating previous elements
        # from $angle_ranges array.
        for my $angle_range ( @{ $angle_ranges } ) {
            $min_angle = $angle_range->[0];
            $max_angle = $angle_range->[1];
            if( $angle >= $min_angle && $angle <= $max_angle ) {
                push @angles, $angle;
                last;
            } elsif( $min_angle == $max_angle ) {
                push @angles, $min_angle;
            }
        }
    }

    return \@angles;
}

#
# Parses query string and generates data structure suited for generating
# rotamers.
# Input:
#     $query_string - query string.
#     E.g. '0..36.0..360.0', '0..18.0..180.0', '0..90.0', 'chi1=0..36.0',
#          'chi1=90.0..90.0, chi2=0.0..10.0..360.0'.
# Output:
#     %angles - data structure describing angles:
#     { 'chi1' => [ 0.0, 1.0, 2.0 ], ... }.
#

sub sample_angles_qs_parsing
{
    my ( $query_string, $in_radians, $small_angle ) = @_;

    $query_string =~ s/\s//g;
    $small_angle = 36.0;

    my %angles;
    for my $angle ( split /,/, $query_string ) {
        my $angle_name;
        my $angle_start;
        my $angle_step;
        my $angle_end;

        if( $angle =~ m/^(\w+)=(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)$/ ) {
            ( $angle_name, $angle_start, $angle_step, $angle_end ) =
                ( $1, $2, $3, $4 );
        } elsif( $angle =~ m/^(\w+)=(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)$/ ) {
            ( $angle_name, $angle_start, $angle_end ) = ( $1, $2, $3 );
        } elsif( $angle =~ m/^(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)$/ ) {
            ( $angle_start, $angle_step, $angle_end ) = ( $1, $2, $3 );
        } elsif( $angle =~ m/^(\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)$/ ) {
            ( $angle_start, $angle_end ) = ( $1, $2 );
        } else {
            die "Syntax '$angle' is incorrect\n"
        }

        $angle_name //= '*';
        $angle_start //= 0.0;
        $angle_step //= $small_angle;
        $angle_end //= 360.0;

        if( $in_radians ) {
            $angles{$angle_name} =
                sample_angles( [ [ $angle_start, $angle_end ] ], $angle_step );
        } else {
            $angles{$angle_name} =
                sample_angles( [ [ $angle_start * $PI / 180.0,
                                   $angle_end * $PI / 180.0 ] ],
                                   $angle_step * $PI / 180.0 );
        }
    }

    return \%angles;
}

1;
