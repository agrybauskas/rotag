package ForceField::Bonded;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( general
                     torsion
                     torsion_new );

use Readonly;

use AtomProperties qw( sort_atom_names );
use Constants qw( $PI );
use ForceField::Parameters;
use Measure qw( dihedral_angle );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Potential functions ----------------------------- #

#
# Calculates bond torsion potential.
#
#                k * ( t_epsilon / 2 ) * ( 1 + cos( n * omega ) )
#
# where:
#     k         - energy weight/adjustment constant;
#     t_epsilon - maximum energy of the peak;
#     n         - number of energy maxima;
#     omega     - dihedral angle.
# Input:
#     $atom_site - atom site data structure;
#     $atom_i_id - target atom id;
#     $parameters - parameters' values.
# Output:
#     $torsion_potential - value of torsion energy potential.
#

sub torsion
{
    my ( $atom_i_id, $parameters ) = @_;

    my ( $t_epsilon, $t_k, $reference_atom_site ) = (
        $parameters->{'t_epsilon'},
        $parameters->{'t_k'},
        $parameters->{'atom_site'}
    );

    $t_epsilon //= 1.0;
    $t_k //= $Parameters::T_K;
    my $t_n //= 3; # FIXME: here the number depends on the hybridization.

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my @connection_ids = @{ $reference_atom_site->{$atom_i_id}{'connections'} };

    my $torsion_potential = 0;
    for my $neighbour_id ( @connection_ids ) {
        my @second_neighbour_ids =
            grep { $atom_i_id ne $_ }
                @{ $reference_atom_site->{$neighbour_id}{'connections'} };

        next if ! @second_neighbour_ids;

        my @second_atom_names =
            @{ sort_atom_names(
                   [ map { $reference_atom_site->{$_}{'label_atom_id'} }
                         @second_neighbour_ids ] ) };
        my ( $second_neighbour_id ) =
            grep { $reference_atom_site->{$_}{'label_atom_id'} eq
                   $second_atom_names[0] }
                 @second_neighbour_ids;

        my @third_neighbour_ids =
            grep { $neighbour_id ne $_ }
                @{ $reference_atom_site->{$second_neighbour_id}{'connections'} };

        next if ! @third_neighbour_ids;

        my @third_atom_names =
            @{ sort_atom_names(
                   [ map { $reference_atom_site->{$_}{'label_atom_id'} }
                         @third_neighbour_ids ] ) };
        my ( $third_neighbour_id ) =
            grep { $reference_atom_site->{$_}{'label_atom_id'} eq
                   $third_atom_names[0] }
                 @third_neighbour_ids;

        my $omega = dihedral_angle(
            [ [ $reference_atom_site->{$third_neighbour_id}{'Cartn_x'},
                $reference_atom_site->{$third_neighbour_id}{'Cartn_y'},
                $reference_atom_site->{$third_neighbour_id}{'Cartn_z'} ],
              [ $reference_atom_site->{$second_neighbour_id}{'Cartn_x'},
                $reference_atom_site->{$second_neighbour_id}{'Cartn_y'},
                $reference_atom_site->{$second_neighbour_id}{'Cartn_z'} ],
              [ $reference_atom_site->{$neighbour_id}{'Cartn_x'},
                $reference_atom_site->{$neighbour_id}{'Cartn_y'},
                $reference_atom_site->{$neighbour_id}{'Cartn_z'} ],
              [ $reference_atom_site->{$atom_i_id}{'Cartn_x'},
                $reference_atom_site->{$atom_i_id}{'Cartn_y'},
                $reference_atom_site->{$atom_i_id}{'Cartn_z'} ] ],
        );

        $torsion_potential +=
            $t_k * ( $t_epsilon / 2 ) * ( 1 + cos( $t_n * $omega ) );
    }

    return $torsion_potential;
}

sub torsion_new
{
    my ( $atom_i_id, $parameters ) = @_;

    my ( $t_k, $reference_atom_site ) = (
        $parameters->{'t_k'},
        $parameters->{'atom_site'}
    );

    $t_k //= $Parameters::T_K;
    my $phase = 3;

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my @connection_ids = @{ $reference_atom_site->{$atom_i_id}{'connections'} };
    my $torsion_potential = 0;
    for my $neighbour_id ( @connection_ids ) {
        my @second_neighbour_ids =
            grep { $atom_i_id ne $_ }
                @{ $reference_atom_site->{$neighbour_id}{'connections'} };

        next if ! @second_neighbour_ids;

        for my $second_neighbour_id ( @second_neighbour_ids ) {
            my @third_neighbour_ids =
                grep { $neighbour_id ne $_ }
                    @{ $reference_atom_site->{$second_neighbour_id}
                                             {'connections'} };

            next if ! @third_neighbour_ids;

            for my $third_neighbour_id ( @third_neighbour_ids ) {
                my $omega = dihedral_angle(
                    [ [ $reference_atom_site->{$third_neighbour_id}{'Cartn_x'},
                        $reference_atom_site->{$third_neighbour_id}{'Cartn_y'},
                        $reference_atom_site->{$third_neighbour_id}{'Cartn_z'} ],
                      [ $reference_atom_site->{$second_neighbour_id}{'Cartn_x'},
                        $reference_atom_site->{$second_neighbour_id}{'Cartn_y'},
                        $reference_atom_site->{$second_neighbour_id}{'Cartn_z'} ],
                      [ $reference_atom_site->{$neighbour_id}{'Cartn_x'},
                        $reference_atom_site->{$neighbour_id}{'Cartn_y'},
                        $reference_atom_site->{$neighbour_id}{'Cartn_z'} ],
                      [ $reference_atom_site->{$atom_i_id}{'Cartn_x'},
                        $reference_atom_site->{$atom_i_id}{'Cartn_y'},
                        $reference_atom_site->{$atom_i_id}{'Cartn_z'} ] ],
                );

                $torsion_potential += $t_k * _torsion( $omega, $phase );
            }
        }
    }

    return $torsion_potential;
}

sub _torsion
{
    my ( $omega, $phase ) = @_;

    if( $omega < ( -$PI / $phase ) || $omega > ( $PI / $phase ) ) {
        return 0
    } else {
        return ( 1 / 2 ) * ( 1 + cos( $phase * $omega ) );
    }
}

sub general
{
    my ( $atom_i, $parameters ) = @_;
    my ( $decompose, $is_optimal ) =
        ( $parameters->{'decompose'}, $parameters->{'is_optimal'} );

    if( $is_optimal ) {
        return 0;
    } else {
        my $torsion = torsion_new( $atom_i->{'id'}, $parameters );

        if( $decompose ) {
            return { 'bonded' => $torsion,
                     'torsion' => $torsion };
        } else {
            return $torsion;
        }
    }
}

1;
