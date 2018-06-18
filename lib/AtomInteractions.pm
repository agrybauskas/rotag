package AtomInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( hard_sphere
                     soft_sphere
                     exponential
                     leonard_jones
                     composite );

use List::Util qw( any max );

use AtomProperties qw( %ATOMS );
use ConnectAtoms qw( connect_atoms
                     grid_box
                     is_connected
                     is_second_neighbour );
use PDBxParser qw( filter );
use LinearAlgebra qw( pi );

# --------------------------- Potential functions ----------------------------- #

#
# Hard sphere potential function. Described as:
#     0,   r_{ij} >= vdw_{i} + vdw_{j}
#     Inf, r_{ij} <  vdw_{i} + vdw_{j}
#
#     where: r - distance between center of atoms;
#            vdw - Van der Waals radius;
#            Inf - infinity;
# Input:
#     $atom_i, $atom_j - atom data structure (see PDBxParser.pm).
# Output:
#     two values: 0 or "Inf" (infinity).
#

sub hard_sphere
{
    my ( $atom_i, $atom_j ) = @_;

    my $vdw_distance =
	$ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
      + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    my $distance =
    	( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
      + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
      + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    if( $distance < $vdw_distance ** 2 ) {
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
#            epsilon - energy coefficient;
#            n - number increases the slope of potential.
# Input:
#     $atom_i, $atom_j - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by soft sphere potential.
#

sub soft_sphere
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my $epsilon =
	$parameters->{'soft_epsilon'} ? $parameters->{'soft_epsilon'} : 1.0;
    my $n =
	$parameters->{'soft_n'} ? $parameters->{'soft_n'} : 12;

    my $vdw_distance =
	$ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
      + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    my $distance =
    	( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
      + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
      + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    if( $distance <= $vdw_distance ** 2 ) {
	return $epsilon * ( $vdw_distance / sqrt( $distance ) )**$n;
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
#            epsilon - energy coefficient;
# Input:
#     $atom_i, $atom_j - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by exponential potential.
#

sub exponential
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my $epsilon =
	$parameters->{'exp_epsilon'} ? $parameters->{'exp_epsilon'} : 1.0;
    my $m =
	$parameters->{'exp_m'} ? $parameters->{'exp_m'} : 1.0;

    my $vdw_distance =
	$ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
      + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    my $distance =
    	( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
      + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
      + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    if( $distance <= $vdw_distance ** 2 ) {
	return $epsilon * exp( - ( sqrt( $distance ) / $vdw_distance ) ** $m );
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
#            epsilon - energy coefficient;
#            sigma - sum of Van der Waals radii of two atoms.
# Input:
#     $atom_i, $atom_j - atom data structure (see PDBxParser.pm).
# Output:
#     value, calculated by Leonard-Jones potential.
#

sub leonard_jones
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my $epsilon =
	$parameters->{'lj_epsilon'} ? $parameters->{'lj_epsilon'} : 1.0;

    my $sigma = # Same as vdw distance.
	$ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
      + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    my $r = # Same as distance.
      sqrt( ( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
          + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
          + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2 );

    return 4 * $epsilon * ( ( $sigma / $r ) ** 12 - ( $sigma / $r ) ** 6 );
}

sub composite
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my $ljones_epsilon =
	$parameters->{'lj_epsilon'} ? $parameters->{'lj_epsilon'} : 1.0;
    my $coulomb_epsilon =
	$parameters->{'c_epsilon'} ? $parameters->{'c_epsilon'} : 0.01;
    my $h_bond_epsilon =
	$parameters->{'h_epsilon'} ? $parameters->{'h_epsilon'} : 1.0;
    my $cutoff_start = # * VdW distance.
	$parameters->{'cutoff_start'} ? $parameters->{'cutoff_start'} : 2.5;
    my $cutoff_end = # * VdW distance.
	$parameters->{'cutoff_end'} ? $parameters->{'cutoff_end'} : 5;

    # Calculates Van der Waals distance of given atoms.
    my $sigma =
    	$ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
      + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    # Calculates squared distance between two atoms.
    my $squared_r =
    	( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
      + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
      + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2;

    # Extracts partial charges.
    my $partial_charge_i =
    	$ATOMS{$atom_i->{'type_symbol'}}{'partial_charge'};
    my $partial_charge_j =
    	$ATOMS{$atom_j->{'type_symbol'}}{'partial_charge'};

    my $hbond; # TODO: hbond will be added, when there will be functions adding
               # hydrogens to atoms.

    if( $squared_r < $cutoff_start * $sigma ) {
	my $leonard_jones =
	    4 * $ljones_epsilon
	      * ( ( ( $sigma / sqrt( $squared_r ) ) ** 12 )
	        - ( ( $sigma / sqrt( $squared_r ) ) ** 6 ) );
	my $coulomb =
	    ( $partial_charge_i * $partial_charge_j ) /
	    ( 4 * $coulomb_epsilon * pi() * sqrt( $squared_r ) );
        return $leonard_jones + $coulomb;
    } elsif( ( $squared_r >= $cutoff_start * $sigma )
	  && ( $squared_r <= $cutoff_end * $sigma ) ) {
	my $leonard_jones =
	    4 * $ljones_epsilon
	      * ( ( ( $sigma / sqrt( $squared_r ) ) ** 12 )
	        - ( ( $sigma / sqrt( $squared_r ) ) ** 6 ) );
	my $coulomb =
	    ( $partial_charge_i * $partial_charge_j ) /
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
