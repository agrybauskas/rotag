package ForceField::Bonded;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( general
                     torsion_components
                     torsion );

use Readonly;

use AtomProperties qw( sort_atom_names );
use Energy;
use ForceField::Parameters;
use Measure qw( dihedral_angle );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------- Potential functions ----------------------------- #

#
# Calculates bond torsion potential.
#
#                t_k * ( epsilon / 2 ) * ( 1 + cos( phase * omega ) )
#
# where:
#     t_k         - energy weight/adjustment constant;
#     epsilon   - maximum energy of the peak;
#     phase     - number of energy maxima;
#     omega     - dihedral angle.
# Input:
#     $parameters - parameters' values.
#     $atom_site - atom site data structure;
#     $atom_i_id - target atom id;
#     $options->{atom_site} - reference atom site.
# Output:
#     $torsion_potential - value of torsion energy potential.
#

sub torsion
{
    my ( $parameters, $atom_i_id, $options ) = @_;

    my ( $reference_atom_site ) = ( $options->{'atom_site'} );

    my $t_k = $parameters->{'_[local]_force_field'}{'t_k'};
    my $torsional = $parameters->{'_[local]_torsional'};

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my $atom_name = $reference_atom_site->{$atom_i_id}{'type_symbol'};
    my @connection_ids = @{ $reference_atom_site->{$atom_i_id}{'connections'} };

    my $torsion_potential_sum = 0;
    my $torsion_potentials =
        torsion_components( $parameters, $atom_i_id, $options );
    for my $torsion_potential ( @{ $torsion_potentials } ) {
        $torsion_potential_sum += $torsion_potential->value;
    }

    return $torsion_potential_sum;
}

sub torsion_components
{
    my ( $parameters, $atom_i_id, $options ) = @_;

    my ( $reference_atom_site ) = ( $options->{'atom_site'} );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $t_k = $parameters->{'_[local]_force_field'}{'t_k'};
    my $torsional = $parameters->{'_[local]_torsional'};

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my $atom_name = $reference_atom_site->{$atom_i_id}{'type_symbol'};
    my @connection_ids = @{ $reference_atom_site->{$atom_i_id}{'connections'} };
    my @torsion_potentials = ();
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
                my $third_atom_name =
                    $reference_atom_site->{$third_neighbour_id}{'type_symbol'};
                my $epsilon = $torsional->{$third_atom_name}{$atom_name}
                                          {'epsilon'};
                my $phase = $torsional->{$third_atom_name}{$atom_name}
                                          {'phase'};
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

                my $torsion_potential;
                if( $omega < ( - $pi / $phase ) || $omega > ( $pi / $phase ) ) {
                    $torsion_potential = 0;
                } else {
                    $torsion_potential =
                        ( $epsilon / 2 ) * ( 1 + cos( $phase * $omega ) );
                }

                my $energy_potential = Energy->new();
                $energy_potential->set_energy(
                    'torsion',
                    [ $atom_i_id, $neighbour_id, $second_neighbour_id,
                      $third_neighbour_id ],
                    $t_k * $torsion_potential
                );
                push @torsion_potentials, $energy_potential;
            }
        }
    }

    return \@torsion_potentials;
}

sub general
{
    my ( $parameters, $atom_i, $options ) = @_;
    my ( $decompose, $is_optimal ) =
        ( $options->{'decompose'}, $options->{'is_optimal'} );

    if( $is_optimal ) {
        return 0;
    } else {
        my $torsion = torsion( $parameters, $atom_i->{'id'}, $options );

        if( $decompose ) {
            return { 'torsion' => $torsion };
        } else {
            return $torsion;
        }
    }
}

1;
