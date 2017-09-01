package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles );

use POSIX;

use lib "./";
use LinearAlgebra qw( pi );

# ---------------------------------- Sampling --------------------------------- #

#
# Generates intervals of values that might be restricted by boundary values.
#

# Produces angle values that are separated by even intervals.
# Input:
#     $angle_ranges - boundary between which angles can be sampled.
#     $small_angle - smallest angle increment.
#     Input has to be in hash form.
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
	2 * pi() / floor( 2 * pi() / $small_angle );
    my @small_angles =
    	map { $_ * $small_angle } 0..( floor( 2 * pi() / $small_angle ) - 1 );

    # Iterates around the circle and adds evenly spaced angles, if they are
    # inside intervals ($angle_ranges).

    for my $angle ( @small_angles ) {
	# TODO: might speed up calculation by eliminating previous elements
	# from $angle_ranges array.
    	for my $angle_range ( @{ $angle_ranges } ) {
    	    $min_angle = $angle_range->[0];
    	    $max_angle = $angle_range->[1];
    	    if( $angle >= $min_angle && $angle <= $max_angle ) {
    		push( @angles, $angle );
		last;
    	    } elsif( $min_angle == $max_angle ) {
		push( @angles, $min_angle );
	    }
    	}
    }

    return \@angles;
}

1;
