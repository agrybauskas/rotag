package ForceField::NonBonded;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( coulomb
                     general
                     h_bond
                     h_bond_explicit
                     h_bond_implicit
                     hard_sphere
                     lennard_jones
                     soft_sphere );

use Carp;
use List::Util qw( any );
use Math::Trig qw( acos
                   asin );
use Readonly;

use ForceField::Parameters;
use Measure qw( bond_angle
                distance_squared );
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
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;

    my ( $r_squared, $sigma ) = ( $options->{'r_squared'}, $options->{'sigma'} );

    my $atom_properties = $parameters->{'_[local]_atom_properties'};

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $atom_properties->{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $atom_properties->{$atom_j->{'type_symbol'}}{'vdw_radius'};

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
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;

    my ( $r_squared, $sigma ) = ( $options->{'r_squared'}, $options->{'sigma'} );

    my $soft_epsilon = $parameters->{'_[local]_force_field'}{'soft_epsilon'};
    my $soft_n = $parameters->{'_[local]_force_field'}{'soft_n'};
    my $atom_properties = $parameters->{'_[local]_atom_properties'};

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $atom_properties->{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $atom_properties->{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared <= $sigma ** 2 ) {
        return $soft_epsilon * ( ( $sigma ** $soft_n ) /
                                 ( $r_squared ** ( $soft_n / 2 ) ) );
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
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;

    my ( $r_squared, $is_optimal  ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
    );

    my $lj_k = $parameters->{'_[local]_force_field'}{'lj_k'};
    my $lennard_jones = $parameters->{'_[local]_lennard_jones'};

    my $lj_epsilon = $lennard_jones->{$atom_i->{'type_symbol'}}
                                     {$atom_j->{'type_symbol'}}
                                     {'epsilon'};

    if( $is_optimal ) {
        return ( -1 ) * $lj_k * $lj_epsilon;
    }

    $r_squared //= distance_squared( $atom_i, $atom_j );

    my $sigma = $lennard_jones->{$atom_i->{'type_symbol'}}
                                {$atom_j->{'type_symbol'}}
                                {'sigma'};

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
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;
    my ( $r_squared, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
    );

    my $c_k = $parameters->{'_[local]_force_field'}{'c_k'};
    my $partial_charge = $parameters->{'_[local]_partial_charge'};

    $r_squared //= distance_squared( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i =
        $partial_charge->{$atom_i->{'label_comp_id'}}{$atom_i->{'label_atom_id'}};
    my $partial_charge_j =
        $partial_charge->{$atom_j->{'label_comp_id'}}{$atom_j->{'label_atom_id'}};

    if( ! defined $partial_charge_i ) {
        confess $atom_i->{'label_atom_id'} . ' atom with id ' . $atom_i->{'id'} .
                ' from ' . $atom_i->{'label_comp_id'} . ' residue does not have '.
                'defined partial charge in force field file' ;
    }
    if( ! defined $partial_charge_j ) {
        confess $atom_j->{'label_atom_id'} . ' atom with id ' . $atom_j->{'id'} .
                ' from ' . $atom_j->{'label_comp_id'} . ' residue does not have '.
                'defined partial charge in force field file' ;
    }

    if( $is_optimal ) {
        if( $partial_charge_i * $partial_charge_j > 0 ) {
            return 0;
        } else {
            # TODO: check if this assumption is true: Lennard-Jones sigma is
            # taken as distance, because Lennard-Jones potential goes faster
            # to infinity than Coulomb.
            my $lennard_jones = $parameters->{'_[local]_lennard_jones'};
            $r_squared = $lennard_jones->{$atom_i->{'type_symbol'}}
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
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;

    my ( $atom_site, $only_implicit ) = (
        $options->{'atom_site'},
        $options->{'only_implicit_h_bond'},
    );

    my $hydrogen_names = $parameters->{'_[local]_hydrogen_names'};

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # BondProperties.
    my @h_bond_heavy_atoms = qw( N O S );
    my @atom_i_hydrogen_names =
        defined $hydrogen_names->{$atom_i->{'label_comp_id'}}
                                 {$atom_i->{'label_atom_id'}} ?
             @{ $hydrogen_names->{$atom_i->{'label_comp_id'}}
                                 {$atom_i->{'label_atom_id'}} } : ();
    my @atom_j_hydrogen_names =
        defined $hydrogen_names->{$atom_j->{'label_comp_id'}}
                                 {$atom_j->{'label_atom_id'}} ?
             @{ $hydrogen_names->{$atom_j->{'label_comp_id'}}
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
            defined $hydrogen_names->{$atom_pair->[0]{'label_comp_id'}}
                                     {$atom_pair->[0]{'label_atom_id'}} ?
                 @{ $hydrogen_names->{$atom_pair->[0]{'label_comp_id'}}
                                     {$atom_pair->[0]{'label_atom_id'}}} : ();

        if( @hydrogen_ids && ! $only_implicit ) {
            for my $hydrogen_id ( @hydrogen_ids ) {
                $h_bond_energy_sum +=
                    h_bond_explicit( $parameters,
                                     $atom_pair->[0],
                                     $atom_site->{$hydrogen_id},
                                     $atom_pair->[1],
                                     $options );
            }
        } elsif( @hydrogen_names ) {
            $h_bond_energy_sum +=
                h_bond_implicit( $parameters,
                                 $atom_pair->[0],
                                 $atom_pair->[1],
                                 $options );
        }
    }

    return $h_bond_energy_sum;
}

#
# Implicit hydrogen bond potential function. Described as explicit, but best
# case scenario is introduced:
#
# h_epsilon * [ 5 * ( sigma / r )**12 - 6 * ( sigma / r )**10 ] * cos( theta )**4
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
#                                   /    /\alpha
#                        (H donor) O ----- O (H acceptor)
#                                          |
#                                          |
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
    my ( $parameters, $donor_atom, $acceptor_atom, $options ) = @_;

    my ( $r_donor_acceptor_squared, $is_optimal, $reference_atom_site ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
        $options->{'atom_site'},
    );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $sp3_angle = $parameters->{'_[local]_constants'}{'sp3_angle'};
    my $sp2_angle = $parameters->{'_[local]_constants'}{'sp2_angle'};
    my $sp_angle = $parameters->{'_[local]_constants'}{'sp_angle'};
    my $h_k = $parameters->{'_[local]_force_field'}->{'h_k'};
    my $atom_properties = $parameters->{'_[local]_atom_properties'};
    my $hydrogen_bond = $parameters->{'_[local]_h_bond'};

    my $r_sigma = $hydrogen_bond->{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon = $hydrogen_bond->{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    if( $is_optimal ) {
        return ( -1 ) * $h_k * $h_epsilon;
    }

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
        $atom_properties->{$donor_atom->{'type_symbol'}}
                          {'covalent_radius'}{'length'}->[$covalent_radius_idx] +
        $atom_properties->{'H'}{'covalent_radius'}{'length'}->[0];
    my $r_acceptor_hydrogen_vdw =
        $atom_properties->{$acceptor_atom->{'type_symbol'}}{'vdw_radius'} +
        $atom_properties->{'H'}{'vdw_radius'};

    # TODO: check for 0, 108.5 angles.
    my $theta;
    if( $r_donor_acceptor_squared < $r_sigma ** 2 ) {
        if( $r_donor_acceptor_squared <
                ( $r_donor_hydrogen + $r_acceptor_hydrogen_vdw ) ** 2 &&
             $r_donor_hydrogen <
                 sqrt( $r_donor_acceptor_squared ) + $r_acceptor_hydrogen_vdw &&
             $r_acceptor_hydrogen_vdw <
                 sqrt( $r_donor_acceptor_squared ) + $r_donor_hydrogen ) {
            $theta = acos(
                ( $r_donor_acceptor_squared -
                  ( $r_donor_hydrogen ** 2 ) -
                  ( $r_acceptor_hydrogen_vdw ** 2 ) ) /
                ( ( -2 ) * $r_donor_hydrogen * $r_acceptor_hydrogen_vdw )
            );
        } else {
            $theta = 0;
        }
    } else {
        # Using both sine and cosine rules we can produce formula that calculates
        # theta of best case scenario:
        #
        #   theta = asin( r_donor_acceptor * sin( alpha ) / r_hydrogen_acceptor )
        #
        my $neighbour_donor_acceptor_angle =
            bond_angle(
                [ [ $reference_atom_site->{$donor_atom->{'connections'}[0]}
                                          {'Cartn_x'},
                    $reference_atom_site->{$donor_atom->{'connections'}[0]}
                                          {'Cartn_y'},
                    $reference_atom_site->{$donor_atom->{'connections'}[0]}
                                          {'Cartn_z'}, ],
                  [ $donor_atom->{'Cartn_x'},
                    $donor_atom->{'Cartn_y'},
                    $donor_atom->{'Cartn_z'}, ],
                  [ $acceptor_atom->{'Cartn_x'},
                    $acceptor_atom->{'Cartn_y'},
                    $acceptor_atom->{'Cartn_z'}, ] ] );

        my $hybridization_angle;
        if( $donor_atom->{'hybridization'} eq 'sp3' ) {
            $hybridization_angle = $sp3_angle;
        } elsif( $donor_atom->{'hybridization'} eq 'sp2' ) {
            $hybridization_angle = $sp2_angle;
        } elsif( $donor_atom->{'hybridization'} eq 'sp' ) {
            $hybridization_angle = $sp_angle;
        }

        my $alpha = abs( $hybridization_angle - $neighbour_donor_acceptor_angle );

        my $r_acceptor_hydrogen_squared =
            $r_donor_hydrogen ** 2 +
            $r_donor_acceptor_squared -
            2 * $r_donor_hydrogen * sqrt( $r_donor_acceptor_squared ) *
            cos( $alpha );

        my $beta = asin(
            sin( $alpha ) * $r_donor_hydrogen /
            sqrt( $r_acceptor_hydrogen_squared )
        );

        $theta = $pi - $alpha - $beta;
    }

    # TODO: study more on what restriction should be on $r_donor_acceptor.
    if( defined $theta && ( $theta >= $pi / 2 ) && ( $theta <=  3 * $pi / 2 ) ) {
        return
            $h_k *
            $h_epsilon *
            ( 5 * ( $r_sigma ** 12 / $r_donor_acceptor_squared ** 6 ) -
              6 * ( $r_sigma ** 10 / $r_donor_acceptor_squared ** 5 ) ) *
            cos( $theta )**4;
    } else {
        return 0;
    }
}

#
# Explicit hydrogen bond potential function. Described as:
#
# h_epsilon * [ 5 * ( sigma / r )**12 - 6 * ( sigma / r )**10 ] * cos( theta )**4
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
    my ( $parameters, $donor_atom, $hydrogen_atom, $acceptor_atom, $options ) = @_;

    my ( $r_donor_acceptor_squared, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'}
    );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $h_k = $parameters->{'_[local]_force_field'}{'h_k'};
    my $hydrogen_bond = $parameters->{'_[local]_h_bond'};

    my $r_sigma = $hydrogen_bond->{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon = $hydrogen_bond->{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    if( $is_optimal ) {
        return (-1) * $h_k * $h_epsilon;
    }

    $r_donor_acceptor_squared //=distance_squared( $donor_atom, $acceptor_atom );
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

    if( ( $theta >= $pi / 2 ) && ( $theta <=  3 * $pi / 2 ) ) {
        return
            $h_k *
            $h_epsilon *
            ( 5 * ( $r_sigma ** 12 / $r_donor_acceptor_squared ** 6 ) -
              6 * ( $r_sigma ** 10 / $r_donor_acceptor_squared ** 5 ) ) *
            cos( $theta )**4;
    } else {
        return 0;
    }
}

#
# Combines Lennard-Jones, Coulomb, hydrogen bond potentials with smoothly
# decreasing cutoff distance.
#
#         general = lennard_jones + coulomb + h_bond * cutoff_function
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

sub general
{
    my ( $parameters, $atom_i, $atom_j, $options ) = @_;

    my ( $r_squared, $sigma, $decompose, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'sigma'},
        $options->{'decompose'}, # Returns hash of the energy function
                                 # component values.
        $options->{'is_optimal'},
    );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $cutoff_start = $parameters->{'_[local]_force_field'}{'cutoff_start'};
    my $cutoff_end = $parameters->{'_[local]_force_field'}{'cutoff_end'};
    my $atom_properties = $parameters->{'_[local]_atom_properties'};

    $decompose //= 0;
    $is_optimal //= 0;

    # Calculates squared distance between two atoms.
    $r_squared //= distance_squared( $atom_i, $atom_j );

    my %options = %{ $options };
    $options{'r_squared'} = $r_squared;

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $atom_properties->{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $atom_properties->{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < ( $cutoff_start * $sigma ) ** 2 ) {
        my $lennard_jones =
            lennard_jones( $parameters, $atom_i, $atom_j,
                           { %options, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $parameters, $atom_i, $atom_j,
                               { %options, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $parameters, $atom_i, $atom_j,
                              { %options, ( 'is_optimal' => $is_optimal ) } );

        if( $decompose ) {
            return { 'lennard_jones' => $lennard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return $lennard_jones + $coulomb + $h_bond;
        }
    } elsif( ( $r_squared >= ( $cutoff_start * $sigma ) ** 2 ) &&
             ( $r_squared <= ( $cutoff_end   * $sigma ) ** 2 ) ) {
        my $lennard_jones =
            lennard_jones( $parameters, $atom_i, $atom_j,
                           { %options, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $parameters, $atom_i, $atom_j,
                               { %options, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $parameters, $atom_i, $atom_j,
                              { %options, ( 'is_optimal' => $is_optimal ) } );
        my $cutoff_function =
            cos( ( $pi * ( sqrt( $r_squared ) - $cutoff_start * $sigma ) ) /
                 ( 2 * ( $cutoff_end * $sigma - $cutoff_start * $sigma ) ) );

        if( $decompose ) {
            return { 'lennard_jones' => $lennard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return ( $lennard_jones + $coulomb + $h_bond ) * $cutoff_function;
        }
    } else {
        if( $decompose ) {
            return { 'lennard_jones' => 0,
                     'coulomb' => 0,
                     'h_bond' => 0 };
        } else {
            return 0;
        }
    }
}

1;
