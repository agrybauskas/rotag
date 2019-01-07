package ForceField::Interactions::Composite;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( general_bonded
                     general_non_bonded);

use ConnectAtoms qw( distance_squared );
use Constants qw( $PI );
use ForceField::Interactions::Bonded qw( torsion );
use ForceField::Interactions::NonBonded qw( coulomb
                                            h_bond
                                            lennard_jones );
use ForceField::Parameters;
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Potential functions ----------------------------- #

#
# Combines Lennard-Jones, Coulomb, hydrogen bond potentials with smoothly
# decreasing cutoff distance.
#
#    general_non_bonded = lennard_jones + coulomb + h_bond * cutoff_function
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

sub general_non_bonded
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

    $cutoff_start //= $Parameters::CUTOFF_START;
    $cutoff_end //= $Parameters::CUTOFF_END;
    $decompose //= 0;
    $is_optimal //= 0;

    # Calculates squared distance between two atoms.
    $r_squared //= distance_squared( $atom_i, $atom_j );

    my %parameters = %{ $parameters };
    $parameters{'r_squared'} = $r_squared;

    # Calculates Van der Waals distance of given atoms.
    $sigma //= $Parameters::ATOMS{$atom_i->{'type_symbol'}}{'vdw_radius'} +
               $Parameters::ATOMS{$atom_j->{'type_symbol'}}{'vdw_radius'};

    if( $r_squared < ( $cutoff_start * $sigma ) ** 2 ) {
        my $lennard_jones =
            lennard_jones( $atom_i, $atom_j,
                           { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $coulomb = coulomb( $atom_i, $atom_j,
                               { %parameters, ( 'is_optimal' => $is_optimal ) } );
        my $h_bond =  h_bond( $atom_i, $atom_j,
                              { %parameters, ( 'is_optimal' => $is_optimal ) } );

        if( $decompose ) {
            return { 'non_bonded' => $lennard_jones + $coulomb + $h_bond,
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
            return { 'non_bonded' =>
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
            return { 'non_bonded' => 0,
                     'lennard_jones' => 0,
                     'coulomb' => 0,
                     'h_bond' => 0 };
        } else {
            return 0;
        }
    }
}

1;
