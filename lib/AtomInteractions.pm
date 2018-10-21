package AtomInteractions;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( composite
                     coulomb
                     h_bond
                     h_bond_explicit
                     h_bond_implicit
                     hard_sphere
                     leonard_jones
                     soft_sphere );

use List::Util qw( any );
use Math::Trig qw( acos );
use Readonly;

use AtomProperties qw( %ATOMS
                       %HYDROGEN_NAMES );
use ConnectAtoms qw( distance
                     distance_squared );
use Constants qw( $CUTOFF_START
                  $CUTOFF_END
                  $COULOMB_K
                  $H_EPSILON
                  $H_SIGMA
                  $LJ_EPSILON
                  $PI
                  $SOFT_EPSILON
                  $SOFT_N );
use Measure qw( bond_angle );
use Version qw( $VERSION );

our $VERSION = $VERSION;

Readonly my $SP3_ANGLE => 109.5 * $PI / 180.0;
Readonly my $SP2_ANGLE => 120.0 * $PI / 180.0;

# --------------------------- Potential functions ----------------------------- #

#
# Hard sphere potential function. Described as:
#                       0,   r_{ij} >= vdw_{i} + vdw_{j}
#                       Inf, r_{ij} <  vdw_{i} + vdw_{j}
#
# where:
#     r   - distance between the center of atoms;
#     vdw - van der Waals radius;
#     Inf - infinity.
#
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     0 or 'Inf' (infinity).
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
        return 'Inf';
    } else {
        return 0;
    }
}

#
# Soft sphere potential function. Described as:
#   epsilon * ( vdw_{i} + vdw_{j} / r_{ij} ) ** n , r_{ij} <= vdw_{i} + vdw_{j}
#   0,                                              r_{ij} >  vdw_{i} + vdw_{j}
#
# where:
#     r       - distance between the centers of the atoms;
#     vdw     - van der Waals radius;
#     epsilon - energy coefficient;
#     n       - number that increases the slope of the potential.
#
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     0 or energy value.
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
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};
    $soft_epsilon //= $SOFT_EPSILON;
    $n //= $SOFT_N;

    if( $r_squared <= $sigma ** 2 ) {
        return $soft_epsilon * ( $sigma / sqrt $r_squared ) ** $n;
    } else {
        return 0;
    }
}

#
# Leonard-Jones potential function. Described as:
#
#          4 * epsilon * [ ( sigma / r ) ** 12 - ( sigma / r ) ** 6 ]
#
# where:
#     r       - distance between center of atoms;
#     epsilon - energy coefficient;
#     sigma   - sum of van der Waals radii of two atoms.
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.
#

sub leonard_jones
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $sigma, $lj_epsilon ) = (
        $parameters->{'r'},
        $parameters->{'sigma'},
        $parameters->{'lj_epsilon'},
    );

    $r //= distance( $atom_i, $atom_j );
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};
    $lj_epsilon //= $LJ_EPSILON;

    return 4 * $lj_epsilon * ( ( $sigma / $r ) ** 12 - ( $sigma / $r ) ** 6 );
}

#
# Coulomb potential function. Described as:
#
#                     coulomb_k * q_{i} * q_{j} / r ** 2
#
# where:
#     r           - distance between center of atoms;
#     coulomb_k   - coulomb coefficient;
#     q_{i,j}     - charge of the particle.
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.
#

sub coulomb
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $coulomb_k ) = ( $parameters->{'r'}, $parameters->{'c_k'} );

    $coulomb_k //= $COULOMB_K;
    $r //= distance( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i = $ATOMS{$atom_i->{'type_symbol'}}{'partial_charge'};
    my $partial_charge_j = $ATOMS{$atom_j->{'type_symbol'}}{'partial_charge'};

    return $coulomb_k * $partial_charge_i * $partial_charge_j / $r ** 2;
}

#
# Calculates hydrogen bond energy by applying implicit and explicit hydrogen
# bond potentials according to the positions and types of atoms.
#
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters->{atom_site} - atom site data structure (see PDBxParser.pm);
#     $parameters->{only_implicit_h_bond} - uses only implicit hydrogen bond
#     potential.
# Output:
#     energy value.

sub h_bond
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $h_epsilon, $atom_site, $only_implicit ) = (
        $parameters->{'h_epsilon'},
        $parameters->{'atom_site'},
        $parameters->{'only_implicit_h_bond'},
    );

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # BondProperties.
    my @h_bond_heavy_atoms = qw( N O F );
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
    if( ! ( ( any { $atom_i->{'type_symbol'} eq $_ } @h_bond_heavy_atoms ) &&
            ( any { $atom_j->{'type_symbol'} eq $_ } @h_bond_heavy_atoms ) ) ||
        ! ( @atom_i_hydrogen_names || @atom_j_hydrogen_names ) ) {
        return 0;
    }

    my $h_bond_energy_sum = 0;

    # Because i and j atoms can be both hydrogen donors and acceptors, two
    # possibilities are explored.
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
                $h_bond_energy_sum +=
                    h_bond_explicit( $atom_pair->[0],
                                     $atom_site->{$hydrogen_id},
                                     $atom_pair->[1],
                                     $parameters );
            }
        } elsif( @hydrogen_names ) {
            $h_bond_energy_sum +=
                h_bond_implicit( $atom_pair->[0], $atom_pair->[1], $parameters );
        }
    }

    return $h_bond_energy_sum;
}

#
# Implicit hydrogen bond potential function. Described as:
#
#                          h_epsilon * cos( theta );
#
# where:
#     h_epsilon - hydrogen bond coefficient;
#     theta     - angle between hydrogen donor, hydrogen and hydrogen acceptor.
#
#                                      H
#                                     /_\theta
#                                    /   \
#                                   /     \
#                        (H donor) O       O (H acceptor)
#
# Input:
#     $donor_atom - donor atom data structure (see PDBxParser.pm);
#     $acceptor atom - acceptor atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.
#

sub h_bond_implicit
{
    my ( $donor_atom, $acceptor_atom, $parameters ) = @_;

    my ( $h_epsilon, $r_sigma ) = (
        $parameters->{'h_epsilon'}, $parameters->{'r_sigma'}
    );

    $r_sigma //= $H_SIGMA;
    $h_epsilon //= $H_EPSILON;

    my $covalent_radius_idx;
    my @hybridizations = qw( sp3 sp2 sp );
    for( my $i = 0; $i <= $#hybridizations; $i++ ) {
        if( $donor_atom->{'hybridization'} eq $hybridizations[$i] ) {
            $covalent_radius_idx = $i;
            last;
        }
    }

    my $r_donor_acceptor = distance( $donor_atom, $acceptor_atom );
    my $r_donor_hydrogen =
        $ATOMS{$donor_atom->{'type_symbol'}}
              {'covalent_radius'}{'length'}->[$covalent_radius_idx] +
        $ATOMS{'H'}{'covalent_radius'}{'length'}->[0];

    # Because there are no information on the position of hydrogen atom, the best
    # case scenario is assumed - the distance of the atoms are optimal. This is
    # done that way in order for explicit function increase energy if needed.
    # That way energy cannot be excluded by the cutoff value.
    my $theta =
        acos( ( $r_donor_acceptor ** 2 - $r_donor_hydrogen ** 2 - $r_sigma ** 2)/
              ( -2 * $r_donor_hydrogen * $r_sigma ) );

    # TODO: study more on what restriction should be on $r_donor_acceptor.
    if( ( $r_donor_acceptor <= $r_donor_hydrogen + $r_sigma ) &&
        ( $theta >= $PI / 2 ) &&
        ( $theta <=  3 * $PI / 2 ) ) {
        return $h_epsilon * cos $theta;
    } else {
        return 0;
    }
}

#
# Explicit hydrogen bond potential function. Described as:
#
# - h_epsilon * [ 5 * ( sigma / r )**12 - 6 * ( sigma / r )**10 ] * cos( theta )
#
# where:
#     r         - distance between center of atoms;
#     h_epsilon - hydrogen bond coefficient;
#     sigma     - optimal distance for hydrogen bond;
#     theta     - angle between hydrogen donor, hydrogen and hydrogen acceptor.
#
#                                      H
#                                     /_\theta
#                                    /   \
#                                   /     \
#                        (H donor) O       O (H acceptor)
#
# Input:
#     $donor_atom - donor atom data structure (see PDBxParser.pm);
#     $hydrogen_atom - hydrogen atom data structure (see PDBxParser.pm);
#     $acceptor_atom - acceptor atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.
#

sub h_bond_explicit
{
    my ( $donor_atom, $hydrogen_atom, $acceptor_atom, $parameters  ) = @_;

    my ( $h_epsilon, $r_sigma ) = (
        $parameters->{'h_epsilon'}, $parameters->{'r_sigma'}
    );
    $r_sigma //= $H_SIGMA;
    $h_epsilon //= $H_EPSILON;

    my $r_acceptor_hydrogen = distance( $hydrogen_atom, $acceptor_atom );
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

    if( ( $theta >= $PI / 2 ) && ( $theta <=  3 * $PI / 2 ) ) {
        return
            ( -1 ) * $h_epsilon *
            ( 5 * ( $r_sigma / $r_acceptor_hydrogen )**12 -
              6 * ( $r_sigma / $r_acceptor_hydrogen )**10 ) *
            cos $theta;
    } else {
        return 0;
    }
}

#
# Combines Leonard-Jones, Coulomb, hydrogen bond potentials with smoothly
# decreasing cutoff distance.
#
#         composite = leonard_jones + coulomb + h_bond * cutoff_function
#
# cutoff_function =
#     cos( ( pi * ( r - cutoff_{start} * sigma ) ) /
#          ( 2 * ( cutoff_{end} * sigma - cutoff_{start} * sigma ) ) )
#
# where:
#     r                  - distance between center of atoms;
#     cutoff_{start,end} - start, end of the cutoff function in amounts of sigma;
#     sigma              - sum of van der Waals radii of two atoms.
#
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.

sub composite
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r, $sigma, $cutoff_start, $cutoff_end, $decompose ) = (
        $parameters->{'r'},
        $parameters->{'sigma'},
        $parameters->{'cutoff_start'}, # * VdW distance.
        $parameters->{'cutoff_end'}, #   * VdW distance.
        $parameters->{'decompose'}, # Returns hash of the energy function
                                    # component values.
    );

    $cutoff_start //= $CUTOFF_START;
    $cutoff_end //= $CUTOFF_END;
    $decompose //= 0;

    # Calculates squared distance between two atoms.
    $r //= sqrt( ( $atom_j->{'Cartn_x'} - $atom_i->{'Cartn_x'} ) ** 2 +
                 ( $atom_j->{'Cartn_y'} - $atom_i->{'Cartn_y'} ) ** 2 +
                 ( $atom_j->{'Cartn_z'} - $atom_i->{'Cartn_z'} ) ** 2 );

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r < $cutoff_start * $sigma ) {
        my $leonard_jones = leonard_jones( $atom_i, $atom_j, $parameters );
        my $coulomb = coulomb( $atom_i, $atom_j, $parameters );
        my $h_bond = h_bond( $atom_i, $atom_j, $parameters );

        if( $decompose ) {
            return { 'composite' => $leonard_jones + $coulomb + $h_bond,
                     'leonard_jones' => $leonard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return $leonard_jones + $coulomb + $h_bond;
        }
    } elsif( ( $r >= $cutoff_start * $sigma ) &&
             ( $r <= $cutoff_end * $sigma ) ) {
        my $leonard_jones = leonard_jones( $atom_i, $atom_j, $parameters );
        my $coulomb = coulomb( $atom_i, $atom_j, $parameters );
        my $h_bond = h_bond( $atom_i, $atom_j, $parameters );
        my $cutoff_function =
            cos( ( $PI * ( $r - $cutoff_start * $sigma ) ) /
                 ( 2 * ( $cutoff_end * $sigma - $cutoff_start * $sigma ) ) );

        if( $decompose ) {
            return { 'composite' =>
                         ( $leonard_jones + $coulomb + $h_bond ) *
                         $cutoff_function,
                     'leonard_jones' => $leonard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return ( $leonard_jones + $coulomb + $h_bond ) * $cutoff_function;
        }
    } else {
        if( $decompose ) {
            return { 'composite' => 0,
                     'leonard_jones' => 0,
                     'coulomb' => 0,
                     'h_bond' => 0 };
        } else {
            return 0;
        }
    }
}

1;
