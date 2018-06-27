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

use AtomProperties qw( %ATOMS
                       %HYDROGEN_NAMES );
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

    my ( $r, $h_epsilon, $atom_site ) = (
        $parameters->{'r'},
        $parameters->{'h_epsilon'},
        $parameters->{'atom_site'}
    );

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # MoleculeProperties.
    my @h_bond_heavy_atoms = ( 'N', 'O', 'F' );
    if( ! ( ( any { $atom_i->{'type_symbol'} eq $_ } @h_bond_heavy_atoms )
         && ( any { $atom_j->{'type_symbol'} eq $_ } @h_bond_heavy_atoms ) ) ) {
        return 0;
    }

    $r //= distance( $atom_i, $atom_j );
    $h_epsilon //= 1.0;

    my $r_i_hydrogen = $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
                     + $ATOMS{'H'}{'vdw_radius'};

    # Calculates the angle (theta) between hydrogen acceptor, hydrogen and
    # hydrogen donor. If there is no information on the position of hydrogens
    # and cannot be determined by hybridization and the quantity of missing
    # hydrogens, the smallest possible is determined by iterating through
    # second neighbours of hydrogen donor and using information about
    # hybridization assign the smallest possible alpha angle.
    #                                    H
    #                                   /_\theta
    #                                  /   \
    #                                 / alpha
    #                      (H donor) O_) _ _ O (H acceptor)
    #
    my $h_bond_energy_sum = 0;

    # TODO: the code looks a bit redundant with i and j cases. Maybe could be
    # refactored.
    my $atom_i_hybridization = $atom_site->{$atom_i->{'id'}}{'hybridization'};
    my $atom_j_hybridization = $atom_site->{$atom_j->{'id'}}{'hybridization'};

    my $atom_i_connection_ids = $atom_site->{$atom_i->{'id'}}{'connections'};
    my $atom_j_connection_ids = $atom_site->{$atom_j->{'id'}}{'connections'};
    my @connection_names_i =
        map { $atom_site->{"$_"}{'label_atom_id'} } @{ $atom_i_connection_ids };
    my @connection_names_j =
        map { $atom_site->{"$_"}{'label_atom_id'} } @{ $atom_j_connection_ids };

    my @atom_i_hydrogen_ids =
        map { $atom_site->{$_}{'type_symbol'} eq 'H' ? $_ : () }
           @{ $atom_i_connection_ids };
    my @atom_j_hydrogen_ids =
        map { $atom_site->{$_}{'type_symbol'} eq 'H' ? $_ : () }
           @{ $atom_j_connection_ids };

    # Checks for missing hydrogens for each atom. If any missing hydrogen is
    # detected, that means clear positions of hydrogens were not predicted,
    # because add_hydrogen() was run with $add_only_clear_positions => 1.
    my $hydrogen_names_i =
        $HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}{$atom_i->{'type_symbol'}};
    my $hydrogen_names_j =
        $HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}{$atom_j->{'type_symbol'}};

    my @missing_hydrogens_i;
    for my $hydrogen_name ( @{ $hydrogen_names_i } ) {
        if( ! grep { /$hydrogen_name/ } @connection_names_i ) {
            push( @missing_hydrogens_i, $hydrogen_name );
        }
    }
    my @missing_hydrogens_j;
    for my $hydrogen_name ( @{ $hydrogen_names_j } ) {
        if( ! grep { /$hydrogen_name/ } @connection_names_j ) {
            push( @missing_hydrogens_j, $hydrogen_name );
        }
    }

    # Because i and j atoms can be both hydrogen donors and acceptors, two
    # possibilities are explored.

    # i-th atom is hydrogen donor and j-th atom - acceptor.

    # j-th atom is hydrogen donor and i-th atom - acceptor.

    # if( ( $theta >= 90 * pi() / 180 ) && ( $theta >= 270 * pi() / 180 ) ) {
    #     return $h_epsilon
    #          * ( 5 * ( $r / $r_i_hydrogen )**12
    #            - 6 * ( $r / $r_i_hydrogen )**10 )
    #          * cos( $theta );
    # } else {
    #     return 0;
    # }

    return $h_bond_energy_sum;
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
