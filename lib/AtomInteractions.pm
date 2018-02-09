package AtomInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( hard_sphere
                     soft_sphere
                     exponential
                     leonard_jones
                     combined );

use List::Util qw( any max );

use AtomProperties qw( %ATOMS );
use ConnectAtoms qw( connect_atoms
                     grid_box
                     is_connected
                     is_second_neighbour );
use PDBxParser qw( filter );
use LinearAlgebra qw( pi );

# ------------------------- Various potential functions ----------------------- #

#
# Hard sphere potential function. Described as:
#     0,   r_{ij} >= vdw_{i} + vdw_{j}
#     Inf, r_{ij} <  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            Inf - infinity;
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     two values: 0 or "Inf" (infinity).
#

sub hard_sphere
{
    my ( $target_atom, $neighbour_atom ) = @_;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance < $vdw_length ** 2 ) {
	return "Inf";
    } else {
	return 0;
    }
}

#
# Soft sphere potential function. Described as:
#     epsilon * ( vdw_{i} + vdw_{j} / r_{ij} ) ** n , r_{ij} <= vdw_{i} + vdw_{j}
#     0, r_{ij} >  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            epsilon - energy coefficient; TODO: should study more about it.
#            n - number increases the slope of potential.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by soft sphere potential.
#

sub soft_sphere
{
    my ( $target_atom, $neighbour_atom, $epsilon, $n ) = @_;

    $epsilon //= 1.0;
    $n //= 12;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance <= $vdw_length ** 2 ) {
	return $epsilon * ( $vdw_length / sqrt( $distance ) )**$n;
    } else {
	return 0;
    }
}


#
# Exponential potential function. Described as:
#     epsilon * exp( - ( r_{ij} / vdw_{i} + vdw_{j} ) ) ** m ,
#        r_{ij} <= vdw_{i} + vdw_{j}
#     0, r_{ij} >  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            epsilon - energy coefficient; TODO: should study more about it.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by exponential potential.
#

sub exponential
{
    my ( $target_atom, $neighbour_atom, $epsilon, $m ) = @_;

    $epsilon //= 1.0;
    $m //= 1.0;

    my $vdw_length =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $distance =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $distance <= $vdw_length ** 2 ) {
	return $epsilon * exp( - ( sqrt( $distance ) / $vdw_length ) ** $m );
    } else {
	return 0;
    }
}

#
# Leonard-Jones potential function. Described as:
#
# 4 * epsilon * [ ( sigma / r ) ** 12 - ( sigma / r ) ** 6 ]
#
#     where: r - distance between center of atoms;
#            epsilon - energy coefficient;  TODO: should study more about it.
#            sigma - sum of Van der Waals radii of two atoms.
# Input:
#     $target_atom, $neighbour_atom - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by Leonard-Jones potential.
#

sub leonard_jones
{
    my ( $target_atom, $neighbour_atom, $epsilon ) = @_;

    $epsilon //= -1.0;

    my $sigma =
	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    my $r =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    if( $r <= ( $sigma ** 2.5 ) ** 2 ) { # TODO:should smoothen function cutoff.
	return 4 * $epsilon * ( ( $sigma / $r ) ** 12 - ( $sigma / $r ) ** 6 );
    } else {
	return 0;
    }
}

sub combined
{
    my ( $target_atom,
	 $neighbour_atom,
	 $ljones_epsilon,
	 $coulomb_epsilon,
	 $h_bond_epsilon,
	 $cutoff_start,
	 $cutoff_end ) = @_;

    # TODO: must check, how combined potential looks in graph.
    $ljones_epsilon //= 1.0;
    $coulomb_epsilon //= 0.01;
    $h_bond_epsilon //= 1.0;
    $cutoff_start //= 2.5; # times VdW distance.
    $cutoff_end //= 5; # times VdW distance.

    # Calculates Van der Waals distance of given atoms.
    my $sigma =
    	$ATOMS{$target_atom->{"type_symbol"}}{"vdw_radius"}
      + $ATOMS{$neighbour_atom->{"type_symbol"}}{"vdw_radius"};

    # Calculates squared distance between two atoms.
    my $squared_r =
    	( $neighbour_atom->{"Cartn_x"} - $target_atom->{"Cartn_x"} ) ** 2
      + ( $neighbour_atom->{"Cartn_y"} - $target_atom->{"Cartn_y"} ) ** 2
      + ( $neighbour_atom->{"Cartn_z"} - $target_atom->{"Cartn_z"} ) ** 2;

    # Extracts partial charges.
    my $target_partial =
    	$ATOMS{$target_atom->{"type_symbol"}}{"partial_charge"};
    my $neighbour_partial =
    	$ATOMS{$neighbour_atom->{"type_symbol"}}{"partial_charge"};

    my $hbond; # TODO: hbond will be added, when there will be functions adding
               # hydrogens to atoms.

    if( $squared_r < $cutoff_start * $sigma ) {
	my $leonard_jones =
	    4 * $ljones_epsilon
	      * ( ( ( $sigma / sqrt( $squared_r ) ) ** 12 )
	        - ( ( $sigma / sqrt( $squared_r ) ) ** 6 ) );
	my $coulomb =
	    ( $target_partial * $neighbour_partial ) /
	    ( 4 * $coulomb_epsilon * pi() * sqrt( $squared_r ) );
        return $leonard_jones + $coulomb;
    } elsif( ( $squared_r >= $cutoff_start * $sigma )
	  && ( $squared_r <= $cutoff_end * $sigma ) ) {
	my $leonard_jones =
	    4 * $ljones_epsilon
	      * ( ( ( $sigma / sqrt( $squared_r ) ) ** 12 )
	        - ( ( $sigma / sqrt( $squared_r ) ) ** 6 ) );
	my $coulomb =
	    ( $target_partial * $neighbour_partial ) /
	    ( 4 * $coulomb_epsilon * pi() * sqrt( $squared_r ) );
	my $cutoff_function =
	    cos( ( pi() * ( sqrt( $squared_r ) - $cutoff_start * $sigma ) ) /
		 ( 2 * ( $cutoff_end * $sigma - $cutoff_start * $sigma ) ) );
	return ( $leonard_jones + $coulomb ) * $cutoff_function;
    } else {
    	return 0;
    }
}

1;
