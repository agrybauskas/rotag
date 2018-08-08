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

use List::Util qw( any );
use Math::Trig qw( acos );

use AtomProperties qw( %ATOMS
                       %HYDROGEN_NAMES );
use ConnectAtoms qw( distance
                     distance_squared );
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

    my ( $r, $coulomb_k ) = ( $parameters->{'r'}, $parameters->{'c_k'} );

    $coulomb_k //= 1.0;
    $r //= distance( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i =
        $ATOMS{$atom_i->{'type_symbol'}}{'partial_charge'};
    my $partial_charge_j =
        $ATOMS{$atom_j->{'type_symbol'}}{'partial_charge'};

    return $coulomb_k * ( $partial_charge_i * $partial_charge_j / $r**2 );
}

sub h_bond
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $h_epsilon, $atom_site, $only_implicit ) = (
        $parameters->{'r'},
        $parameters->{'h_epsilon'},
        $parameters->{'atom_site'},
        $parameters->{'only_implicit_h_bond'}
    );

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # MoleculeProperties.
    # TODO: read about the situations when there are hydrogen atom in the
    # middle of two hydrogen acceptors.
    my @h_bond_heavy_atoms = ( 'N', 'O', 'F' );
    my @atom_i_hydrogen_names =
        defined $HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                               {$atom_i->{'label_atom_id'}} ?
             @{ $HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                               {$atom_i->{'label_atom_id'}} } : ();
    my @atom_j_hydrogen_names =
        defined $HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
                               {$atom_j->{'label_atom_id'}} ?
             @{ $HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
                               {$atom_j->{'label_atom_id'}} } : ();

    # Exits early if there are no hydrogen bond donor-acceptor combinations.
    if( ! ( ( any { $atom_i->{'type_symbol'} eq $_ } @h_bond_heavy_atoms )
         && ( any { $atom_j->{'type_symbol'} eq $_ } @h_bond_heavy_atoms ) )
     || ! ( @atom_i_hydrogen_names || @atom_j_hydrogen_names ) ) {
        return 0;
    }

    $r //= distance( $atom_i, $atom_j );
    $h_epsilon //= -2.5e+05;

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

    # Because i and j atoms can be both hydrogen donors and acceptors, two
    # possibilities are explored.
    my @h_bonds; # List of hashes are used: ( { 'theta' => <float>,  } )

    for my $atom_pair ( [ $atom_i, $atom_j ], [ $atom_j, $atom_i ] ) {
        my @hydrogen_ids =
            map { $atom_site->{$_}{'type_symbol'} eq 'H' ? $_ : () }
               @{ $atom_pair->[0]{'connections'} };
        my @hydrogen_names =
            defined $HYDROGEN_NAMES{$atom_pair->[0]{'label_comp_id'}}
                                   {$atom_pair->[0]{'label_atom_id'}} ?
                 @{ $HYDROGEN_NAMES{$atom_pair->[0]{'label_comp_id'}}
                                   {$atom_pair->[0]{'label_atom_id'}}} : ();

        if( @hydrogen_ids && ! $only_implicit ) {
            for my $hydrogen_id ( @hydrogen_ids ) {
                push( @h_bonds,
                      h_bond_explicit( $atom_i,
                                       $atom_site->{$hydrogen_id},
                                       $atom_j ) );
            }
        } elsif( @hydrogen_names ) {
            push( @h_bonds,
                  h_bond_implicit( $atom_site, $atom_i, $atom_j ) );
        }
    }

    # Calculates the sum of all hydrogen bonds.
    for my $h_bond ( @h_bonds ) {
        if( ( $h_bond->{'theta'} >=   90 * pi() / 180 )
         && ( $h_bond->{'theta'} <=  270 * pi() / 180 ) ) {
            $h_bond_energy_sum +=
                $h_epsilon * ( 5 * ( $h_bond->{'r_donor_hydrogen'} / $r )**12
                             - 6 * ( $h_bond->{'r_donor_hydrogen'} / $r )**10 )
              * cos( $h_bond->{'theta'} );
        }
    }

    return $h_bond_energy_sum;
}

sub h_bond_implicit
{
    my ( $atom_site, $donor_atom, $acceptor_atom ) = @_;

    my $covalent_radius_idx;
    my @hybridizations = ( 'sp3', 'sp2', 'sp' );
    for( my $i = 0; $i <= $#hybridizations; $i++ ) {
        if( $donor_atom->{'hybridization'} eq $hybridizations[$i] ) {
            $covalent_radius_idx = $i;
            last;
        }
    }

    my $r_donor_hydrogen =
        $ATOMS{$donor_atom->{'type_symbol'}}
              {'covalent_radius'}{'length'}->[$covalent_radius_idx]
      + $ATOMS{'H'}{'covalent_radius'}{'length'}->[0];

    # Determines smallest and most favourable angle which theta will be
    # calculated from. The smaller alpha angle, the greater theta is.
    my $alpha;
    if( $donor_atom->{'hybridization'} eq 'sp3' ) {
        $alpha = 109.5 * pi() / 180; # TODO: consider electron pairs?
    } elsif( $donor_atom->{'hybridization'} eq 'sp2' ) {
        $alpha = 120 * pi() / 180;
    }

    my $alpha_delta;
    for my $donor_connection_id ( @{ $donor_atom->{'connections'} } ) {
        my $alpha_delta_local =
            bond_angle(
                [ [ $acceptor_atom->{'Cartn_x'},
                    $acceptor_atom->{'Cartn_y'},
                    $acceptor_atom->{'Cartn_z'} ],
                  [ $donor_atom->{'Cartn_x'},
                    $donor_atom->{'Cartn_y'},
                    $donor_atom->{'Cartn_z'} ],
                  [ $atom_site->{$donor_connection_id}{'Cartn_x'},
                    $atom_site->{$donor_connection_id}{'Cartn_y'},
                    $atom_site->{$donor_connection_id}{'Cartn_z'} ] ] );

        if( ! defined $alpha_delta ) {
            $alpha_delta = $alpha_delta_local
        } elsif( $alpha_delta_local < $alpha_delta ) {
            $alpha_delta = $alpha_delta_local;
        }
    }

    $alpha = $alpha - $alpha_delta;

    my $r = distance( $donor_atom, $acceptor_atom );
    my $theta = acos(
        ( $r_donor_hydrogen - $r * cos( $alpha ) )
      / sqrt( $r_donor_hydrogen**2
            + $r**2
            - 2 * $r_donor_hydrogen * $r * cos( $alpha ) )
    );

    return {'theta' => $theta, 'r_donor_hydrogen' => $r_donor_hydrogen};
}

sub h_bond_explicit
{
    my ( $donor_atom, $hydrogen_atom, $acceptor_atom  ) = @_;

    my $theta = bond_angle(
        [ [ $donor_atom->{'Cartn_x'},
            $donor_atom->{'Cartn_y'},
            $donor_atom->{'Cartn_z'} ],
          [ $hydrogen_atom->{'Cartn_x'},
            $hydrogen_atom->{'Cartn_y'},
            $hydrogen_atom->{'Cartn_z'} ],
          [ $acceptor_atom->{'Cartn_x'},
            $acceptor_atom->{'Cartn_y'},
            $acceptor_atom->{'Cartn_z'} ] ] );
    my $r_donor_hydrogen = distance( $donor_atom, $hydrogen_atom );

    return { 'theta' => $theta, 'r_donor_hydrogen' => $r_donor_hydrogen };
}

sub composite
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $sigma, $cutoff_start, $cutoff_end ) = (
        $parameters->{'r'},
        $parameters->{'sigma'},
        $parameters->{'cutoff_start'}, # * VdW distance.
        $parameters->{'cutoff_end'}, #   * VdW distance.
    );

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
