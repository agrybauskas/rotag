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
use Math::Trig;

use AtomProperties qw( %ATOMS
                       %HYDROGEN_NAMES );
use ConnectAtoms qw( connect_atoms
                     distance
                     distance_squared
                     grid_box
                     is_connected
                     is_second_neighbour );
use PDBxParser qw( filter );
# use LinearAlgebra qw( pi );
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
    # TODO: read about the situations when there are hydrogen atom in the
    # middle of two hydrogen acceptors.
    my @h_bond_heavy_atoms = ( 'N', 'O', 'F' );
    my @atom_i_hydrogen_names =
        defined $HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                               {$atom_i->{'type_symbol'}} ?
                $HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                               {$atom_i->{'type_symbol'}} : ();
    my @atom_j_hydrogen_names =
        defined $HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
                               {$atom_j->{'type_symbol'}} ?
                $HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
                               {$atom_j->{'type_symbol'}} : ();

    # Exits early if there are no hydrogen bond donor-acceptor combinations.
    if( ! ( ( any { $atom_i->{'type_symbol'} eq $_ } @h_bond_heavy_atoms )
         && ( any { $atom_j->{'type_symbol'} eq $_ } @h_bond_heavy_atoms ) )
     || ! ( @atom_i_hydrogen_names || @atom_j_hydrogen_names ) ) {
        return 0;
    }

    $r //= distance( $atom_i, $atom_j );
    $h_epsilon //= 1.0;

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

    # Because i and j atoms can be both hydrogen donors and acceptors, two
    # possibilities are explored.
    my @atom_i_hydrogen_ids =
        map { $atom_site->{$_}{'type_symbol'} eq 'H' ? $_ : () }
           @{ $atom_i_connection_ids };
    my @atom_j_hydrogen_ids =
        map { $atom_site->{$_}{'type_symbol'} eq 'H' ? $_ : () }
           @{ $atom_j_connection_ids };

    my @atom_i_coord =
        ( $atom_i->{'Cartn_x'}, $atom_i->{'Cartn_y'}, $atom_i->{'Cartn_z'} );
    my @atom_j_coord =
        ( $atom_j->{'Cartn_x'}, $atom_j->{'Cartn_y'}, $atom_j->{'Cartn_z'} );

    my @h_bonds; # List of hashes are used: ( { 'theta' => <float>,  } )

    # i-th atom is hydrogen donor and j-th atom - acceptor.
    if( @atom_i_hydrogen_ids ) {
        for my $atom_i_hydrogen_id ( @atom_i_hydrogen_ids ) {
            my $theta = bond_angle(
                [ \@atom_i_coord,
                  [ $atom_site->{$atom_i_hydrogen_id}{'Cartn_x'},
                    $atom_site->{$atom_i_hydrogen_id}{'Cartn_y'},
                    $atom_site->{$atom_i_hydrogen_id}{'Cartn_z'} ],
                  \@atom_j_coord ] );
            my $r_i_hydrogen =
                distance( $atom_i, $atom_site->{$atom_i_hydrogen_id} );
            push( @h_bonds,
                  { 'theta' => $theta, 'r_i_hydrogen' => $r_i_hydrogen } );
        }
    } elsif( @atom_i_hydrogen_names ) { # Hydrogens are generalized here. See
        my $r_i_hydrogen =              # comments above.
            $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
          + $ATOMS{'H'}{'vdw_radius'};

        # Determines smallest and most favourable angle which theta will be
        # calculated from. The smaller alpha angle, the greater theta is.
        my $alpha;
        if( $atom_i_hybridization eq 'sp3' ) {
            $alpha = 109.5 * pi() / 180; # TODO: should consider electron pairs?
        } elsif( $atom_i_hybridization eq 'sp2' ) {
            $alpha = 120 * pi() / 180; # TODO: should be hardcoded?
        } # sp hybridization is absent, because you can determine clear geometry.

        # TODO: check if it is possible to be in this loop if there are hydrogens
        # connected?
        my $alpha_delta = 0; # Value that should be substracted from $alpha.
        for my $atom_i_connection_id ( @{ $atom_i_connection_ids } ) {
            my $alpha_delta_local =
                bond_angle(
                    [ \@atom_i_coord,
                      \@atom_j_coord,
                      [ $atom_site->{$atom_i_connection_id}{'Cartn_x'},
                        $atom_site->{$atom_i_connection_id}{'Cartn_y'},
                        $atom_site->{$atom_i_connection_id}{'Cartn_z'} ] ] );

            # TODO: check if this statement below is actually true. For now,
            # it seems like it.
            if( $alpha_delta_local < $alpha_delta ) {
                $alpha_delta = $alpha_delta_local;
            }
        }
        $alpha = $alpha - $alpha_delta;

        my $theta = acos(
            ( $r_i_hydrogen - $r * cos( $alpha ) )
          / sqrt( $r_i_hydrogen**2 + $r**2 - 2 * $r_i_hydrogen * $r * cos( $alpha ) )
        );

        push( @h_bonds, { 'theta' => $theta, 'r_i_hydrogen' => $r_i_hydrogen } );
    }

    # j-th atom is hydrogen donor and i-th atom - acceptor.

    for my $h_bond ( @h_bonds ) {
        if( ( $h_bond->{'theta'} >= 90 * pi() / 180 )
         && ( $h_bond->{'theta'} >= 270 * pi() / 180 ) ) {
            $h_bond_energy_sum =
                $h_bond_energy_sum
              + $h_epsilon * ( 5 * ( $r / $h_bond->{'r_i_hydrogen'} )**12
                             - 6 * ( $r / $h_bond->{'r_i_hydrogen'} )**10 )
                * cos( $h_bond->{'theta'} );
        }
    }

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
