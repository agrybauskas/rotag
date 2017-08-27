package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles );

use lib "./";
use LinearAlgebra qw( epsilon );

# ---------------------------------- Sampling --------------------------------- #

#
# Generates intervals of values that might be restricted by boundary values.
#

# Produces angle values that are separated by even intervals.
# Input:
#     $angle_ranges - boundary between which angles can be sampled.
#     $angle_df - smallest angle increment.
#     Input has to be in hash form.
# Output:
#     @angles - sampled angles.
#
#

sub sample_angles
{
    my ( $angle_ranges, $angle_df ) = @_;

    my @angles;
    my $min_angle;
    my $max_angle;
    my $center_angle;
    my $current_angle;

    for my $angle_range ( @{ $angle_ranges } ) {
	# Because $angle_range might be smaller than $angle_df or is not
	# devisable by $angle_df, additional angles have to be added outside
	# the range so, boundary conditions could be tested.

	# Finds center of the range and starts incrementing angles by $angle_df.
	$min_angle = $angle_range->[0];
	$max_angle = $angle_range->[1];
	$center_angle = ( $min_angle + $max_angle ) / 2;
	push( @angles, $center_angle ); # Center point is added to array of
	                                # angles.

	# Increases angles  positively.
	$current_angle = $center_angle;
	while( $current_angle <= $max_angle ) {
	    $current_angle += $angle_df;
	    push( @angles, $current_angle );
	}

	# Increases angles  negatively.
	$current_angle = $center_angle - $angle_df;
	while( $current_angle >= $min_angle ) {
	    $current_angle -= $angle_df;
	    push( @angles, $current_angle );
	}
    }

    return \@angles;
}

1;
