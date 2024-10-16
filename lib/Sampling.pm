package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles
                     sample_angles_qs_parsing );

use POSIX;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Sampling ---------------------------------- #

#
# Produces angle values that are separated by even intervals.
# Input:
#     $angle_ranges - boundary between which angles can be sampled.
#     $small_angle - smallest angle increment.
# Output:
#     @angles - sampled angles.
#
#

sub sample_angles
{
    my ( $parameters, $angle_ranges, $small_angle, $angle_phase_shift, $rand_count ) = @_;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    $angle_phase_shift //= - $pi;

    my @angles;
    my $min_angle;
    my $max_angle;

    if( defined $rand_count ) {
        # TODO: add random sampling.
    } else {
        # Devides full circle (2*pi) into even intervals by $small_angle value.
        $small_angle = # Adjusts angle so, it could be devided evenly.
            2 * $pi / floor( 2 * $pi / $small_angle );
        my @small_angles =
            map { $_ * $small_angle + $angle_phase_shift }
                ( 0..( floor( 2 * $pi / $small_angle ) - 1 ) );

        # Iterates around the circle and adds evenly spaced angles, if they are
        # inside intervals ($angle_ranges).
        for my $angle ( @small_angles ) {
            # TODO: might speed up calculation by eliminating previous elements
            # from $angle_ranges array.
            for my $angle_range ( @{ $angle_ranges } ) {
                $min_angle = $angle_range->[0];
                $max_angle = $angle_range->[1];
                if( $angle >= $min_angle && $angle <= $max_angle ) {
                    push @angles, $angle;
                    last;
                } elsif( $min_angle == $max_angle ) {
                    push @angles, $min_angle;
                }
            }
        }
    }

    return \@angles;
}

#
# Parses query strings and generates data structure suited for generating
# rotamers.
# Input:
#     $query_strings - query strings.
#     E.g. '0..36.0..360.0', '0..18.0..180.0', '0..90.0', 'chi1=0..36.0',
#          'chi1=90.0..90.0, chi2=0.0..10.0..360.0'.
# Output:
#     %angles - data structure describing angles:
#     { 'chi1' => [ 0.0, 1.0, 2.0 ], ... }.
#

sub sample_angles_qs_parsing
{
    my ( $parameters, $query_strings, $in_radians, $small_angle ) = @_;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $dihedral_angle_restraints =
        $parameters->{'_[local]_dihedral_angle_restraints'};
    my $rotatable_residue_names =
        $parameters->{'_[local]_rotatable_residue_names'};

    $query_strings =~ s/\s//g;
    $small_angle = 36.0;

    my %angles;
    for my $residue_name ( sort keys %{ $dihedral_angle_restraints } ) {
        for my $angle_name ( sort keys %{ $dihedral_angle_restraints->{$residue_name} } ) {
            if( $residue_name eq '.' ) {
                $residue_name = '*'
            }
            if( $angle_name eq '.' ) {
                $angle_name = '*'
            }
            my ( $angle_start, $angle_step, $angle_end ) =
                retrieve_dihedral_angle_params( $dihedral_angle_restraints,
                                                $residue_name,
                                                $angle_name,
                                                [ 'range_from', 'step', 'range_to' ] );

            if( $in_radians ) {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( $parameters, [ [ $angle_start, $angle_end ] ],
                                   $angle_step );
            } else {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( $parameters,
                                   [ [ $angle_start * $pi / 180.0,
                                       $angle_end * $pi / 180.0 ] ],
                                   $angle_step * $pi / 180.0 );
            }
        }
    }

    # Query overwrites on top.
    if( $query_strings ) {
        undef %angles;
    }
    for my $query_string ( split /;/, $query_strings ) {
        # HACK: it should be moved to grammar module and generalized.
        my $residue_names;
        my $angle_string;

        my $residue_names_regexp = join '|', @{ $rotatable_residue_names };
        my @query_string_decomposed = split /:/, $query_string;
        if( scalar @query_string_decomposed == 2 ) {
            if( $query_string =~ m/^((?:${residue_names_regexp})(?:,(?:${residue_names_regexp}))*):(.+)$/i ) {
                $residue_names = [ split /,/, uc( $1 ) ];
                $angle_string = $2;
            } else {
                die "Syntax '$query_string' is incorrect\n"
            }
        } elsif( scalar @query_string_decomposed == 1 ) {
            $angle_string = $query_string;
        } else {
            die "Syntax '$query_string' is incorrect\n"
        }

        $residue_names //= [ "*" ];
        $angle_string //= "";

        for my $angle ( split /,/, $angle_string ) {
            my $angle_name;
            my $angle_start;
            my $angle_step;
            my $angle_end;

            if( $angle =~ m/^(\w+)=(-?\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_name, $angle_start, $angle_step, $angle_end ) =
                    ( $1, $2, $3, $4 );
            } elsif( $angle =~ m/^(\w+)=(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_name, $angle_start, $angle_end ) = ( $1, $2, $3 );
            } elsif( $angle =~ m/^(\w+)=(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_name, $angle_step ) = ( $1, $2 );
            } elsif( $angle =~ m/^(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_step ) = ( $1 );
            } elsif( $angle =~ m/^(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_start, $angle_step, $angle_end ) = ( $1, $2, $3 );
            } elsif( $angle =~ m/^(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $angle_start, $angle_end ) = ( $1, $2 );
            } else {
                die "Syntax '$angle' is incorrect\n"
            }

            $angle_name //= '*';
            $angle_start //= - 180.0;
            $angle_step //= $small_angle;
            $angle_end //= 180.0;

            for my $residue_name ( @{ $residue_names } ) {
                if( $in_radians ) {
                    $angles{$residue_name}{$angle_name} =
                        sample_angles( $parameters,
                                       [ [ $angle_start, $angle_end ] ],
                                       $angle_step );
                } else {
                    $angles{$residue_name}{$angle_name} =
                        sample_angles( $parameters,
                                       [ [ $angle_start * $pi / 180.0,
                                           $angle_end * $pi / 180.0 ] ],
                                       $angle_step * $pi / 180.0 );
                }
            }
        }
    }

    return \%angles;
}

sub retrieve_dihedral_angle_params
{
    my ( $dihedral_angle_restraints, $residue_name, $angle_name, $params ) = @_;

    my %params = ();
    for my $param ( @{ $params } ) {
        my $angle_specific =
            $dihedral_angle_restraints->{$residue_name}{$angle_name}{$param};
        my $residue_specific =
            $dihedral_angle_restraints->{$residue_name}{'.'}{$param};
        my $nonspecific =
            $dihedral_angle_restraints->{'.'}{'.'}{$param};

        if( defined $angle_specific && $angle_specific ne '.' ) {
            $params{$param} = $angle_specific;
        } elsif( defined $residue_specific && $residue_specific ne '.' ) {
            $params{$param} = $residue_specific;
        } else {
            $params{$param} = $nonspecific;
        }
    }

    return map { $params{$_} } @{ $params };
}

1;
