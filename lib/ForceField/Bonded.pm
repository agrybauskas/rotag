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

    my @connection_ids =
        exists $reference_atom_site->{$atom_i_id}{'connections'} ?
        ( @{ $reference_atom_site->{$atom_i_id}{'connections'} } ) : ();

    return 0 if ! @connection_ids;

    my $t_k = $parameters->{'_[local]_force_field'}{'t_k'};
    my $torsional = $parameters->{'_[local]_torsional'};

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my $atom_name = $reference_atom_site->{$atom_i_id}{'type_symbol'};

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

    my ( $reference_atom_site, $debug ) =
        ( $options->{'atom_site'}, $options->{'debug'} );

    $debug //= 0;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $t_k = $parameters->{'_[local]_force_field'}{'t_k'};
    my $torsional = $parameters->{'_[local]_torsional'};
    my $torsional_atom_names = $parameters->{'_[local]_torsional_atom_names'};

    # Determines all dihedral angles by searching third neighbours following the
    # connections.
    my $atom_name = $reference_atom_site->{$atom_i_id}{'label_atom_id'};
    my $residue_name = $reference_atom_site->{$atom_i_id}{'label_comp_id'};
    my @connection_ids =
        exists $reference_atom_site->{$atom_i_id}{'connections'} ?
        ( @{ $reference_atom_site->{$atom_i_id}{'connections'} } ) : ();
    my @torsion_potentials = ();
    for my $neighbour_id ( @connection_ids ) {
        my $neighbour_atom_name =
            $reference_atom_site->{$neighbour_id}{'label_atom_id'};
        my @second_neighbour_ids =
            grep { $atom_i_id ne $_ }
                @{ $reference_atom_site->{$neighbour_id}{'connections'} };

        next if ! @second_neighbour_ids;

        for my $second_neighbour_id ( @second_neighbour_ids ) {
            my $second_atom_name =
                $reference_atom_site->{$second_neighbour_id}{'label_atom_id'};
            my @third_neighbour_ids =
                grep { $neighbour_id ne $_ }
                    @{ $reference_atom_site->{$second_neighbour_id}
                                             {'connections'} };

            next if ! @third_neighbour_ids;

            for my $third_neighbour_id ( @third_neighbour_ids ) {
                my $alt_atom_name =
                    $torsional_atom_names->{$residue_name}
                                           {$atom_name}[0];
                my $alt_neighbour_atom_name =
                    $torsional_atom_names->{$residue_name}
                                           {$neighbour_atom_name}[0];
                my $alt_second_atom_name =
                    $torsional_atom_names->{$residue_name}
                                           {$second_atom_name}[0];
                my $alt_third_atom_name =
                    $torsional_atom_names->{$residue_name}
                                           {$reference_atom_site->
                                                {$third_neighbour_id}
                                                {'label_atom_id'}}[0];

                my $epsilon;
                my $phase;
                my $gamma;

                if( defined $alt_atom_name &&
                    defined $alt_neighbour_atom_name &&
                    defined $alt_second_atom_name &&
                    defined $alt_third_atom_name ) {
                    my @torsion_angle_keys = (
                        "$alt_atom_name,$alt_neighbour_atom_name," .
                        "$alt_second_atom_name,$alt_third_atom_name",
                        "$alt_third_atom_name,$alt_second_atom_name," .
                        "$alt_neighbour_atom_name,$alt_atom_name",
                        "?,$alt_neighbour_atom_name,$alt_second_atom_name,?",
                        "?,$alt_second_atom_name,$alt_neighbour_atom_name,?",
                    );
                    for my $torsion_angle_key ( @torsion_angle_keys ) {
                        if( defined $torsional->{$torsion_angle_key} ) {
                            if( $debug ) {
                                print STDERR
                                    "$atom_i_id,$neighbour_id,".
                                    "$second_neighbour_id,$third_neighbour_id," .
                                    "$atom_name,$neighbour_atom_name," .
                                    "$second_atom_name," .
                                    $reference_atom_site->{$third_neighbour_id}
                                                          {'label_atom_id'} .",".
                                    "$torsion_angle_key,";
                            }

                            ( $epsilon, $phase, $gamma ) = (
                                $torsional->{$torsion_angle_key}{'epsilon'},
                                $torsional->{$torsion_angle_key}{'phase'},
                                $pi * $torsional->{$torsion_angle_key}{'gamma'} / 180.0,
                            );

                            last;
                        }
                    }
                }

                if( ! defined $epsilon && ! defined $phase &&
                    ! defined $gamma ) {
                    ( $epsilon, $phase, $gamma ) = ( 0.0, 1.0, 0.0 );
                }

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

                if( $debug ) {
                    print STDERR sprintf( '%.3f', 180.0 * $omega / $pi ), "\n";
                }

                my $torsion_potential;
                if( $omega < - ( $pi - $gamma ) / $phase ||
                    $omega >   ( $pi + $gamma ) / $phase ) {
                    $torsion_potential = 0;
                } else {
                    $torsion_potential =
                        ( $epsilon / 2 ) * ( 1 + cos( $phase * $omega - $gamma ) );
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
