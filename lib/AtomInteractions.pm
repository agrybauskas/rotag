package AtomInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( hard_sphere
                     soft_sphere
                     leonard_jones
                     coulomb
                     h_bond
                     composite );

use List::Util qw( any max );

use AtomProperties qw( %ATOMS );
use ConnectAtoms qw( connect_atoms
                     distance
                     distance_squared
                     grid_box
                     is_connected
                     is_second_neighbour );
use PDBxParser qw( filter );
use LinearAlgebra qw( pi );
use Measure qw( bond_angle );

# --------------------------- Potential functions ----------------------------- #

# TODO: rewrite function descriptions.

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
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r_squared, $sigma, $soft_epsilon, $n ) = (
        $parameters->{'r_squared'},
        $parameters->{'sigma'},
    );

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < $sigma ** 2 ) {
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

    my ( $r_squared, $sigma, $soft_epsilon, $n ) = (
        $parameters->{'r_squared'},
        $parameters->{'sigma'},
        $parameters->{'soft_epsilon'},
        $parameters->{'soft_n'},
    );

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};
    $soft_epsilon //= 1.0;
    $n //= 12;

    if( $r_squared <= $sigma ** 2 ) {
	return $soft_epsilon * ( $sigma / sqrt( $r_squared ) )**$n;
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

    my ( $r, $sigma, $lj_epsilon ) = (
        $parameters->{'r'},
        $parameters->{'sigma'},
        $parameters->{'lj_epsilon'}
    );

    $r //= distance( $atom_i, $atom_j );
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};
    $lj_epsilon //= 1.0;

    return 4 * $lj_epsilon * ( ( $sigma / $r ) ** 12 - ( $sigma / $r ) ** 6 );
}

sub coulomb
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $coulomb_epsilon ) =
        ( $parameters->{'r'}, $parameters->{'c_epsilon'} );

    $coulomb_epsilon //= 0.1;
    $r //= distance( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i =
    	$ATOMS{$atom_i->{'type_symbol'}}{'partial_charge'};
    my $partial_charge_j =
    	$ATOMS{$atom_j->{'type_symbol'}}{'partial_charge'};

    return ( $partial_charge_i * $partial_charge_j ) /
	   ( 4 * $coulomb_epsilon * pi() * $r );
}

sub h_bond
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $h_epsilon, $connections ) = (
        $parameters->{'r'},
        $parameters->{'h_epsilon'},
        $parameters->{'connections'},
    );

    $r //= distance( $atom_i, $atom_j );
    $h_epsilon //= 1.0;
    my $r_i_hydrogen = $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
                     + $ATOMS{'H'}{'vdw_radius'};

    # my $theta = bond_angle(
    #     [ $atom_i->{'Cartn_x'}, $atom_i->{'Cartn_y'}, $atom_i->{'Cartn_z'} ],
    #     [ $atom_j->{'Cartn_x'}, $atom_j->{'Cartn_y'}, $atom_j->{'Cartn_z'} ],
    #     [ $atom_h->{'Cartn_x'}, $atom_h->{'Cartn_y'}, $atom_h->{'Cartn_z'} ]
    # );

    # if( ( $theta >= 90 * pi() / 180 ) && ( $theta >= 270 * pi() / 180 ) ) {
    #     return $h_epsilon
    #          * ( ( $r / $r_i_hydrogen )**12 - ( $r / $r_i_hydrogen )**10 );
    # } else {
    #     return 0;
    # }

    return 0;
}

sub composite
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $sigma, $lj_epsilon, $coulomb_epsilon,
         $h_bond_epsilon, $cutoff_start, $cutoff_end ) = (
        $parameters->{'r'},
        $parameters->{'sigma'},
        $parameters->{'lj_epsilon'},
        $parameters->{'c_epsilon'},
        $parameters->{'h_epsilon'},
        $parameters->{'cutoff_start'}, # * VdW distance.
        $parameters->{'cutoff_end'}, #   * VdW distance.
    );

    $lj_epsilon //= 1.0;
    $coulomb_epsilon //= 0.1;
    $h_bond_epsilon //= 1.0;
    $cutoff_start //= 2.5;
    $cutoff_end //= 5.0;

    # Calculates squared distance between two atoms.
    $r //= sqrt( ( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2
                 + ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2
                 + ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2 );

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    # For the sake of not repeating r and sigma calculations, they should be
    # passed as the parameters.
    $parameters->{'r'} = $r;
    $parameters->{'sigma'} = $sigma;

    if( $r < $cutoff_start * $sigma ) {
        my $leonard_jones = leonard_jones( $atom_i, $atom_j, $parameters );
        my $coulomb = coulomb( $atom_i, $atom_j, $parameters );
        my $h_bond = h_bond( $atom_i, $atom_j, $parameters );
        return $leonard_jones + $coulomb + $h_bond;
    } elsif( ( $r >= $cutoff_start * $sigma )
          && ( $r <= $cutoff_end * $sigma ) ) {
        my $leonard_jones = leonard_jones( $atom_i, $atom_j, $parameters );
        my $coulomb = coulomb( $atom_i, $atom_j, $parameters );
        my $h_bond = h_bond( $atom_i, $atom_j, $parameters );
        my $cutoff_function =
            cos( ( pi() * ( $r - $cutoff_start * $sigma ) ) /
        	 ( 2 * ( $cutoff_end * $sigma - $cutoff_start * $sigma ) ) );
        return ( $leonard_jones + $coulomb + $h_bond ) * $cutoff_function;
    } else {
    	return 0;
    }
}

1;
