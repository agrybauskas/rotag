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

use ConnectAtoms qw( distance_squared );
use ForceField::Parameters;
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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;

    my ( $r_squared, $sigma ) = ( $options->{'r_squared'}, $options->{'sigma'} );

    my $ATOM_PROPERTIES = $PARAMETERS->{'_[local]_atom_properties'};

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $ATOM_PROPERTIES->{$atom_i->{'type_symbol'}}{'vdw_radius'}
             + $ATOM_PROPERTIES->{$atom_j->{'type_symbol'}}{'vdw_radius'};

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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;

    my ( $r_squared, $sigma ) = ( $options->{'r_squared'}, $options->{'sigma'} );

    my $SOFT_EPSILON = $PARAMETERS->{'_[local]_force_field'}{'soft_epsilon'};
    my $SOFT_N = $PARAMETERS->{'_[local]_force_field'}{'soft_n'};
    my $ATOM_PROPERTIES = $PARAMETERS->{'_[local]_atom_properties'};

    $r_squared //= distance_squared( $atom_i, $atom_j );
    $sigma //= $ATOM_PROPERTIES->{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $ATOM_PROPERTIES->{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared <= $sigma ** 2 ) {
        return $SOFT_EPSILON * ( ( $sigma ** $SOFT_N ) /
                                 ( $r_squared ** ( $SOFT_N / 2 ) ) );
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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;

    my ( $r_squared, $is_optimal  ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
    );

    my $LJ_K = $PARAMETERS->{'_[local]_force_field'}{'lj_k'};
    my $LENNARD_JONES = $PARAMETERS->{'_[local]_lennard_jones'};

    my $lj_epsilon = $LENNARD_JONES->{$atom_i->{'type_symbol'}}
                                     {$atom_j->{'type_symbol'}}
                                     {'epsilon'};

    if( $is_optimal ) {
        return (-1) * $LJ_K * $lj_epsilon;
    }

    $r_squared //= distance_squared( $atom_i, $atom_j );

    my $sigma = $LENNARD_JONES->{$atom_i->{'type_symbol'}}
                                {$atom_j->{'type_symbol'}}
                                {'sigma'};

    return 4 * $LJ_K * $lj_epsilon * ( ( $sigma ** 12 / $r_squared ** 6 ) -
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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;
    my ( $r_squared, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
    );

    my $C_K = $PARAMETERS->{'_[local]_force_field'}{'c_k'};
    my $PARTIAL_CHARGE = $PARAMETERS->{'_[local]_partial_charge'};

    $r_squared //= distance_squared( $atom_i, $atom_j );

    # Extracts partial charges.
    my $partial_charge_i =
        $PARTIAL_CHARGE->{$atom_i->{'label_comp_id'}}{$atom_i->{'label_atom_id'}};
    my $partial_charge_j =
        $PARTIAL_CHARGE->{$atom_j->{'label_comp_id'}}{$atom_j->{'label_atom_id'}};

    if( ! defined $partial_charge_i ) {
        confess $atom_i->{'label_atom_id'} . 'atom with id ' . $atom_i->{'id'} .
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
            my $LENNARD_JONES = $PARAMETERS->{'_[local]_lennard_jones'};
            $r_squared = $LENNARD_JONES->{$atom_i->{'type_symbol'}}
                                         {$atom_j->{'type_symbol'}}
                                         {'sigma'} ** 2;
        }
    }

    return $C_K * $partial_charge_i * $partial_charge_j / $r_squared;
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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;

    my ( $atom_site, $only_implicit ) = (
        $options->{'atom_site'},
        $options->{'only_implicit_h_bond'},
    );

    my $HYDROGEN_NAMES = $PARAMETERS->{'_[local]_hydrogen_names'};

    # TODO: should not be hardcoded - maybe stored in AtomProperties or
    # BondProperties.
    my @h_bond_heavy_atoms = qw( N O S );
    my @atom_i_hydrogen_names =
        defined $HYDROGEN_NAMES->{$atom_i->{'label_comp_id'}}
                                 {$atom_i->{'label_atom_id'}} ?
             @{ $HYDROGEN_NAMES->{$atom_i->{'label_comp_id'}}
                                 {$atom_i->{'label_atom_id'}} } : ();
    my @atom_j_hydrogen_names =
        defined $HYDROGEN_NAMES->{$atom_j->{'label_comp_id'}}
                                 {$atom_j->{'label_atom_id'}} ?
             @{ $HYDROGEN_NAMES->{$atom_j->{'label_comp_id'}}
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
            defined $HYDROGEN_NAMES->{$atom_pair->[0]{'label_comp_id'}}
                                     {$atom_pair->[0]{'label_atom_id'}} ?
                 @{ $HYDROGEN_NAMES->{$atom_pair->[0]{'label_comp_id'}}
                                     {$atom_pair->[0]{'label_atom_id'}}} : ();

        if( @hydrogen_ids && ! $only_implicit ) {
            for my $hydrogen_id ( @hydrogen_ids ) {
                $h_bond_energy_sum +=
                    h_bond_explicit( $atom_pair->[0],
                                     $atom_site->{$hydrogen_id},
                                     $atom_pair->[1],
                                     $PARAMETERS,
                                     $options );
            }
        } elsif( @hydrogen_names ) {
            $h_bond_energy_sum +=
                h_bond_implicit( $atom_pair->[0],
                                 $atom_pair->[1],
                                 $PARAMETERS,
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
    my ( $donor_atom, $acceptor_atom, $PARAMETERS, $options ) = @_;

    my ( $r_donor_acceptor_squared, $is_optimal, $reference_atom_site ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'},
        $options->{'atom_site'},
    );

    my $PI = $PARAMETERS->{'_[local]_constants'}{'pi'};
    my $SP3_ANGLE = $PARAMETERS->{'_[local]_constants'}{'sp3_angle'};
    my $SP2_ANGLE = $PARAMETERS->{'_[local]_constants'}{'sp2_angle'};
    my $SP_ANGLE = $PARAMETERS->{'_[local]_constants'}{'sp_angle'};
    my $H_K = $PARAMETERS->{'_[local]_force_field'}->{'h_k'};
    my $ATOM_PROPERTIES = $PARAMETERS->{'_[local]_atom_properties'};
    my $HYDROGEN_BOND = $PARAMETERS->{'_[local]_h_bond'};

    my $r_sigma = $HYDROGEN_BOND->{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon = $HYDROGEN_BOND->{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    if( $is_optimal ) {
        return (-1) * $H_K * $h_epsilon;
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
        $ATOM_PROPERTIES->{$donor_atom->{'type_symbol'}}
                          {'covalent_radius'}{'length'}->[$covalent_radius_idx] +
        $ATOM_PROPERTIES->{'H'}{'covalent_radius'}{'length'}->[0];
    my $r_acceptor_hydrogen_vdw =
        $ATOM_PROPERTIES->{$acceptor_atom->{'type_symbol'}}{'vdw_radius'} +
        $ATOM_PROPERTIES->{'H'}{'vdw_radius'};

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
            $hybridization_angle = $SP3_ANGLE;
        } elsif( $donor_atom->{'hybridization'} eq 'sp2' ) {
            $hybridization_angle = $SP2_ANGLE;
        } elsif( $donor_atom->{'hybridization'} eq 'sp' ) {
            $hybridization_angle = $SP_ANGLE;
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

        $theta = $PI - $alpha - $beta;
    }

    # TODO: study more on what restriction should be on $r_donor_acceptor.
    if( defined $theta && ( $theta >= $PI / 2 ) && ( $theta <=  3 * $PI / 2 ) ) {
        return
            $H_K *
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
    my ( $donor_atom, $hydrogen_atom, $acceptor_atom, $PARAMETERS, $options  ) = @_;

    my ( $r_donor_acceptor_squared, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'is_optimal'}
    );

    my $PI = $PARAMETERS->{'_[local]_constants'}{'pi'};
    my $H_K = $PARAMETERS->{'_[local]_force_field'}{'h_k'};
    my $HYDROGEN_BOND = $PARAMETERS->{'_[local]_h_bond'};

    my $r_sigma = $HYDROGEN_BOND->{$acceptor_atom->{'type_symbol'}}{'sigma'};
    my $h_epsilon = $HYDROGEN_BOND->{$acceptor_atom->{'type_symbol'}}{'epsilon'};

    if( $is_optimal ) {
        return (-1) * $H_K * $h_epsilon;
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

    if( ( $theta >= $PI / 2 ) && ( $theta <=  3 * $PI / 2 ) ) {
        return
            $H_K *
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
    my ( $atom_i, $atom_j, $PARAMETERS, $options ) = @_;

    my ( $r_squared, $sigma, $decompose, $is_optimal ) = (
        $options->{'r_squared'},
        $options->{'sigma'},
        $options->{'decompose'}, # Returns hash of the energy function
                                 # component values.
        $options->{'is_optimal'},
    );

    my $PI = $PARAMETERS->{'_[local]_constants'}{'pi'};
    my $CUTOFF_START = $PARAMETERS->{'_[local]_force_field'}{'cutoff_start'};
    my $CUTOFF_END = $PARAMETERS->{'_[local]_force_field'}{'cutoff_end'};
    my $ATOM_PROPERTIES = $PARAMETERS->{'_[local]_atom_properties'};

    $decompose //= 0;
    $is_optimal //= 0;

    # Calculates squared distance between two atoms.
    $r_squared //= distance_squared( $atom_i, $atom_j );

    my %options = %{ $options };
    $options{'r_squared'} = $r_squared;

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $ATOM_PROPERTIES->{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $ATOM_PROPERTIES->{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < ( $CUTOFF_START * $sigma ) ** 2 ) {
        my $lennard_jones =
            lennard_jones( $atom_i, $atom_j, $PARAMETERS,
                           { %options, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $atom_i, $atom_j, $PARAMETERS,
                               { %options, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $atom_i, $atom_j, $PARAMETERS,
                              { %options, ( 'is_optimal' => $is_optimal ) } );

        if( $decompose ) {
            return { 'lennard_jones' => $lennard_jones,
                     'coulomb' => $coulomb,
                     'h_bond' => $h_bond };
        } else {
            return $lennard_jones + $coulomb + $h_bond;
        }
    } elsif( ( $r_squared >= ( $CUTOFF_START * $sigma ) ** 2 ) &&
             ( $r_squared <= ( $CUTOFF_END   * $sigma ) ** 2 ) ) {
        my $lennard_jones =
            lennard_jones( $atom_i, $atom_j, $PARAMETERS,
                           { %options, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $atom_i, $atom_j, $PARAMETERS,
                               { %options, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $atom_i, $atom_j, $PARAMETERS,
                              { %options, ( 'is_optimal' => $is_optimal ) } );
        my $cutoff_function =
            cos( ( $PI * ( sqrt( $r_squared ) - $CUTOFF_START * $sigma ) ) /
                 ( 2 * ( $CUTOFF_END * $sigma - $CUTOFF_START * $sigma ) ) );

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
