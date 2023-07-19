package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles
                     sample_bond_parameters
                     sample_angles_qs_parsing );

use POSIX;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Sampling ---------------------------------- #

#
# Produces angle values that are separated by even intervals.
# Input:
#     $parameters - data structure from Parameters.pm;
#     $angle_ranges - boundary between which angles can be sampled;
#     $angle_count - sampling count.
# Output:
#     sampled angles.
#
#

sub sample_angles
{
    my ( $parameters, $angle_ranges, $angle_count ) = @_;
    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    return sample_bond_parameters( $angle_ranges, $angle_count, - $pi, 1, 0 );
}

#
# Produces bond parameter values that are separated by even intervals.
# Input:
#     $bond_parameter_ranges - boundary between which bond parameters can be
#     sampled;
#     $sampling_count - sampling count;
#     $bond_parameter_shift - used when the values have to be moved by certain
#     value;
#     $inclusive_start - inclusive range start and related to mathematical
#     notation (0, 2];
#     $inclusive_end -
#     $inclusive_end - inclusive range end and related to mathematical notation
#     [0, 2).
# Output:
#     \@bond_parameter_values - sampled bond parameter values.
#

sub sample_bond_parameters
{
    my ( $bond_parameter_ranges, $sampling_count, $bond_parameter_shift,
         $inclusive_start, $inclusive_end ) = @_;

    return [] if ! $sampling_count;

    $bond_parameter_shift //= 0;
    $inclusive_start //= 1;
    $inclusive_end //= 1;

    my @bond_parameter_values;

    my $min_value = $bond_parameter_ranges->[0][0];
    my $max_value = $bond_parameter_ranges->[0][1];

    my $updated_sampling_count =
        $sampling_count +
        ( $inclusive_start && $inclusive_end && $sampling_count > 1 ? -1 : 0 );
    my $small_change =
        ( $max_value - $min_value ) / $updated_sampling_count;

    my $sampling_adjustment = 0;
    if( ! $inclusive_start && ! $inclusive_end ) {
        $sampling_adjustment = $small_change / 2;
    } elsif( ! $inclusive_start ) {
        $sampling_adjustment = $small_change;
    }

    @bond_parameter_values =
        map { $min_value + $_ * $small_change + $sampling_adjustment +
              $bond_parameter_shift }
            ( 0..$sampling_count - 1 );

    return \@bond_parameter_values;
}

#
# Parses query string and generates data structure suited for generating
# rotamers.
# Input:
#     $query_string - query string.
#     E.g. '0..36.0..360.0', '0..18.0..180.0', '0..90.0', 'chi1=0..36.0',
#          'chi1=90.0..90.0, chi2=0.0..10.0..360.0'.
# Output:
#     %angles - data structure describing angles:
#     { 'chi1' => [ 0.0, 1.0, 2.0 ], ... }.
#

sub sample_angles_qs_parsing
{
    my ( $parameters, $query_string, $in_radians, $small_angle ) = @_;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $dihedral_angle_restraints =
        $parameters->{'_[local]_dihedral_angle_restraints'};

    $query_string =~ s/\s//g;
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

            my $angle_count = int( ( $angle_end - $angle_start ) / $angle_step );

            if( $in_radians ) {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( $parameters, [ [ $angle_start, $angle_end ] ],
                                   $angle_count );
            } else {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( $parameters,
                                   [ [ $angle_start * $pi / 180.0,
                                       $angle_end * $pi / 180.0 ] ],
                                   $angle_count );
            }
        }
    }

    # Query overwrites on top.
    if( $query_string ) {
        undef %angles;
    }
    for my $angle ( split /,/, $query_string ) {
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
        }else {
            die "Syntax '$angle' is incorrect\n"
        }

        $angle_name //= '*';
        $angle_start //= - 180.0;
        $angle_step //= $small_angle;
        $angle_end //= 180.0;

        my $angle_count = int( ( $angle_end - $angle_start ) / $angle_step );

        if( $in_radians ) {
            $angles{'*'}{$angle_name} =
                sample_angles( $parameters, [ [ $angle_start, $angle_end ] ],
                               $angle_count );
        } else {
            $angles{'*'}{$angle_name} =
                sample_angles( $parameters,
                               [ [ $angle_start * $pi / 180.0,
                                   $angle_end * $pi / 180.0 ] ],
                               $angle_count );
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
