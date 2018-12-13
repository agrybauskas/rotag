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
                     lennard_jones
                     soft_sphere );

use List::Util qw( any );
use Math::Trig qw( acos );
use Readonly;

use ConnectAtoms qw( distance_squared );
use Constants qw( $PI );
use ForceField::General;
use Measure qw( bond_angle );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Potential functions ----------------------------- #

#
# Hard sphere potential function. Described as:
#                       0,   r_ij >= vdw_i + vdw_j
#                       Inf, r_ij <  vdw_i + vdw_j
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
    $sigma //= $General::ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $General::ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < $sigma ** 2 ) {
        return 'Inf';
    } else {
        return 0;
    }
}

#
# Soft sphere potential function. Described as:
#   epsilon * ( vdw_i + vdw_j / r_ij ) ^ n , r_ij <= vdw_i + vdw_j
#   0,                                       r_ij >  vdw_i + vdw_j
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
    $sigma //= $General::ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $General::ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};
    $soft_epsilon //= $General::SOFT_EPSILON;
    $n //= $General::SOFT_N;

    if( $r_squared <= $sigma ** 2 ) {
        return $soft_epsilon * ( ( $sigma ** $n ) /
                                 ( $r_squared ** ( $n / 2 ) ) );
    } else {
        return 0;
    }
}

#
# Lennard-Jones potential function. Described as:
#
#          4 * k * epsilon * [ ( sigma / r ) ** 12 - ( sigma / r ) ** 6 ]
#
# where:
#     r       - distance between center of atoms;
#     epsilon - energy coefficient;
#     sigma   - sum of van der Waals radii of two atoms;
#     k       - weight/adjustment constant.
# Input:
#     $atom_{i,j} - atom data structure (see PDBxParser.pm);
#     $parameters - hash of parameters' values.
# Output:
#     energy value.
#

sub lennard_jones
{
    my ( $atom_i, $atom_j, $parameters ) = @_;

    my ( $r_squared, $lj_k, $is_optimal  ) = (
        $parameters->{'r_squared'},
        $parameters->{'lj_k'},
        $parameters->{'is_optimal'},
    );

    $lj_k //= $General::LJ_K;

    if( $is_optimal ) {
        return (-1) * $lj_k * $General::LENNARD_JONES{$atom_i->{'type_symbol'}}
                                                     {$atom_j->{'type_symbol'}}
                                                     {'epsilon'};
    }

    $r_squared //= distance_squared( $atom_i, $atom_j );

    my $sigma = $General::LENNARD_JONES{$atom_i->{'type_symbol'}}
                                       {$atom_j->{'type_symbol'}}
                                       {'sigma'};
    my $lj_epsilon = $General::LENNARD_JONES{$atom_i->{'type_symbol'}}
                                            {$atom_j->{'type_symbol'}}
                                            {'epsilon'};

    return 4 * $lj_k * $lj_epsilon * ( ( $sigma ** 12 / $r_squared ** 6 ) -
                                       ( $sigma ** 6  / $r_squared ** 3 ) );
}

#
# Coulomb potential function. Described as:
#
#                         k * q_i * q_j / r ** 2
#
# where:
#     r           - distance between center of atoms;
#     k           - weight/adjustment constant;
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
    my ( $r_squared, $c_k, $is_optimal ) = (
        $parameters->{'r_squared'},
        $parameters->{'c_k'},
        $parameters->{'is_optimal'},
    );

    $c_k //= $General::C_K;
    $r_squared //= distance_squared( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i =
        $General::PARTIAL_CHARGE{$atom_i->{'label_comp_id'}}
                                {$atom_i->{'label_atom_id'}};
    my $partial_charge_j =
        $General::PARTIAL_CHARGE{$atom_j->{'label_comp_id'}}
                                {$atom_j->{'label_atom_id'}};

    if( $is_optimal ) {
        if( $partial_charge_i * $partial_charge_j > 0 ) {
            return 0;
        } else {
            # TODO: check if this assumption is true: Lennard-Jones sigma is
            # taken as distance, because Lennard-Jones potential goes faster
            # to infinity than Coulomb.
            $r_squared = $General::LENNARD_JONES{$atom_i->{'type_symbol'}}
                                                {$atom_j->{'type_symbol'}}
                                                {'sigma'} ** 2;
        }
    }

    return $c_k * $partial_charge_i * $partial_charge_j / $r_squared;
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

    my ( $atom_site, $only_implicit ) = (
        $parameters->{'atom_site'},
        $parameters->{'only_implicit_h_bond'},
    );

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # BondProperties.
    my @h_bond_heavy_atoms = qw( N O S );
    my @atom_i_hydrogen_names =
        defined $General::HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                                        {$atom_i->{'label_atom_id'}} ?
             @{ $General::HYDROGEN_NAMES{$atom_i->{'label_comp_id'}}
                                        {$atom_i->{'label_atom_id'}} } : ();
    my @atom_j_hydrogen_names =
        defined $General::HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
                                        {$atom_j->{'label_atom_id'}} ?
             @{ $General::HYDROGEN_NAMES{$atom_j->{'label_comp_id'}}
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
        if( ! exists $atom_pair->[0]{'connections'} ||
            ! defined $atom_pair->[0]{'connections'} ||
            ! @{ $atom_pair->[0]{'connections'} } ) {
            next;
        };

        my @hydrogen_ids =
            map { ( exists $atom_site->{$_} &&
                    $atom_site->{$_}{'type_symbol'} eq 'H' ) ? $_ : () }
               @{ $atom_pair->[0]{'connections'} };
        my @hydrogen_names =
            defined $General::HYDROGEN_NAMES{$atom_pair->[0]{'label_comp_id'}}
                                            {$atom_pair->[0]{'label_atom_id'}} ?
                 @{ $General::HYDROGEN_NAMES{$atom_pair->[0]{'label_comp_id'}}
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

    my ( $r_donor_acceptor_squared, $h_k, $is_optimal ) = (
        $parameters->{'r_squared'},
        $parameters->{'h_k'},
        $parameters->{'is_optimal'},
    );

    $h_k //= $General::H_K;

    my $r_sigma =
        $General::HYDROGEN_BOND{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon =
        $General::HYDROGEN_BOND{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    my $covalent_radius_idx;
    my @hybridizations = qw( sp3 sp2 sp );
    for( my $i = 0; $i <= $#hybridizations; $i++ ) {
        if( $donor_atom->{'hybridization'} eq $hybridizations[$i] ) {
            $covalent_radius_idx = $i;
            last;
        }
    }

    $r_donor_acceptor_squared //= distance_squared( $donor_atom, $acceptor_atom);
    my $r_donor_hydrogen =
        $General::ATOMS{$donor_atom->{'type_symbol'}}
                       {'covalent_radius'}{'length'}->[$covalent_radius_idx] +
        $General::ATOMS{'H'}{'covalent_radius'}{'length'}->[0];

    # Because there are no information on the position of hydrogen atom, the best
    # case scenario is assumed - the distance of the atoms are optimal. This is
    # done that way in order for explicit function increase energy if needed.
    # That way energy cannot be excluded by the cutoff value.
    my $theta =
        acos( ( $r_donor_acceptor_squared - $r_donor_hydrogen**2 - $r_sigma**2)/
              ( -2 * $r_donor_hydrogen * $r_sigma ) );

    # TODO: study more on what restriction should be on $r_donor_acceptor.
    if( ( $r_donor_acceptor_squared <= ( $r_donor_hydrogen + $r_sigma ) ** 2 ) &&
        ( $theta >= $PI / 2 ) &&
        ( $theta <=  3 * $PI / 2 ) ) {
        return $h_k * $h_epsilon * cos $theta;
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

    my ( $h_k, $is_optimal ) = (
        $parameters->{'h_k'},
        $parameters->{'is_optimal'},
    );

    $h_k //= $General::H_K;

    my $r_sigma =
        $General::HYDROGEN_BOND{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon =
        $General::HYDROGEN_BOND{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    my $r_acceptor_hydrogen_squared =
        distance_squared( $hydrogen_atom, $acceptor_atom );
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
            ( -1 ) *
            $h_k *
            $h_epsilon *
            ( 5 * ( $r_sigma ** 12 / $r_acceptor_hydrogen_squared ** 6 ) -
              6 * ( $r_sigma ** 10 / $r_acceptor_hydrogen_squared ** 5 ) ) *
            cos $theta;
    } else {
        return 0;
    }
}

#
# Combines Lennard-Jones, Coulomb, hydrogen bond potentials with smoothly
# decreasing cutoff distance.
#
#         composite = lennard_jones + coulomb + h_bond * cutoff_function
#
# cutoff_function =
#     cos( ( pi * ( r - cutoff_start * sigma ) ) /
#          ( 2  * ( cutoff_end * sigma - cutoff_start * sigma ) ) )
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

    my ( $r_squared, $sigma, $cutoff_start, $cutoff_end, $decompose,
         $is_optimal ) = (
        $parameters->{'r_squared'},
        $parameters->{'sigma'},
        $parameters->{'cutoff_start'}, # * VdW distance.
        $parameters->{'cutoff_end'}, # * VdW distance.
        $parameters->{'decompose'}, # Returns hash of the energy function
                                    # component values.
        $parameters->{'is_optimal'},
    );

    $cutoff_start //= $General::CUTOFF_START;
    $cutoff_end //= $General::CUTOFF_END;
    $decompose //= 0;
    $is_optimal //= 0;

    # Calculates squared distance between two atoms.
    $r_squared //= distance_squared( $atom_i, $atom_j );

    my %parameters = %{ $parameters };
    $parameters{'r_squared'} = $r_squared;

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $General::ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $General::ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < ( $cutoff_start * $sigma ) ** 2 ) {
        my $lennard_jones =
            lennard_jones( $atom_i, $atom_j,
                           { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $atom_i, $atom_j,
                               { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $atom_i, $atom_j,
                              { %parameters, ( 'is_optimal' => $is_optimal ) } );

        if( $decompose ) {
            return { 'composite' => $lennard_jones + $coulomb + $h_bond,
                     'lennard_jones' => $lennard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return $lennard_jones + $coulomb + $h_bond;
        }
    } elsif( ( $r_squared >= ( $cutoff_start * $sigma ) ** 2 ) &&
             ( $r_squared <= ( $cutoff_end   * $sigma ) ** 2 ) ) {
        my $lennard_jones =
            lennard_jones( $atom_i, $atom_j,
                           { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $atom_i, $atom_j,
                               { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $atom_i, $atom_j,
                              { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $cutoff_function =
            cos( ( $PI * ( sqrt( $r_squared ) - $cutoff_start * $sigma ) ) /
                 ( 2 * ( $cutoff_end * $sigma - $cutoff_start * $sigma ) ) );

        if( $decompose ) {
            return { 'composite' =>
                         ( $lennard_jones + $coulomb + $h_bond ) *
                         $cutoff_function,
                     'lennard_jones' => $lennard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return ( $lennard_jones + $coulomb + $h_bond ) * $cutoff_function;
        }
    } else {
        if( $decompose ) {
            return { 'composite' => 0,
                     'lennard_jones' => 0,
                     'coulomb' => 0,
                     'h_bond' => 0 };
        } else {
            return 0;
        }
    }
}

1;
