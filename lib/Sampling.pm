package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles
                     sample_bond_parameters
                     sample_bond_parameters_qs_parsing );

use POSIX;

use BondParameters qw( detect_bond_parameter_type );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Sampling ---------------------------------- #

#
# Produces angle values that are separated by even intervals.
# Input:
#     $angle_ranges - boundary between which angles can be sampled;
#     $angle_count - sampling count.
# Output:
#     sampled angles.
#
#

sub sample_angles
{
    my ( $angle_ranges, $angle_count ) = @_;
    return sample_bond_parameters( $angle_ranges, $angle_count, 1, 0 );
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
    my ( $bond_parameter_ranges, $sampling_count, $inclusive_start,
         $inclusive_end ) = @_;

    return [] if ! $sampling_count;

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
        map { $min_value + $_ * $small_change + $sampling_adjustment }
            ( 0..$sampling_count - 1 );

    return \@bond_parameter_values;
}

#
# Parses query strings and generates data structure suited for generating
# rotamers.
# Input:
#     $query_strings - query strings.
#     E.g. '0..36.0..360.0', '0..18.0..180.0', '0..90.0', 'chi1=0..36.0',
#          'chi1=90.0..90.0, chi2=0.0..10.0..360.0'.
# Output:
#     %bond_parameters - data structure describing bond_parameters:
#     { 'SER' => { 'chi1' => [ 0.0, 1.0, 2.0 ], ... } }.
#

sub sample_bond_parameters_qs_parsing
{
    my ( $parameters, $query_strings, $options ) = @_;

    my ( $in_radians, $legacy_grammar ) =
        ( $options->{'in_radians'}, $options->{'legacy_grammar'} );

    $legacy_grammar //= 0;

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $dihedral_angle_restraints =
        $parameters->{'_[local]_bond_parameter_restraints'};
    my $rotatable_residue_names =
        $parameters->{'_[local]_rotatable_residue_names'};

    $query_strings =~ s/\s//g;

    # Dihedral angle calculations are first preditermined and then updated. For
    # bond angle and bond length parameters, there are no defaults and they
    # have to be declared as it is not desirable to calculate each bond angle
    # and bond length changes.
    my %bond_parameters;
    for my $residue_name ( sort keys %{ $dihedral_angle_restraints } ) {
        for my $dihedral_angle_name (
            sort keys %{ $dihedral_angle_restraints->{$residue_name} } ) {
            my ( $dihedral_angle_start,
                 $dihedral_angle_step,
                 $dihedral_angle_end ) =
                retrieve_dihedral_angle_params( $dihedral_angle_restraints,
                                                $residue_name,
                                                $dihedral_angle_name,
                                                [ 'range_from',
                                                  'step',
                                                  'range_to' ] );

            my $dihedral_angle_count =
                int( ( $dihedral_angle_end - $dihedral_angle_start ) /
                     $dihedral_angle_step );

            if( $in_radians ) {
                $bond_parameters{$residue_name}{$dihedral_angle_name} = {
                    'values' =>
                        sample_angles( [ [ $dihedral_angle_start,
                                           $dihedral_angle_end ] ],
                                       $dihedral_angle_count ),
                    'type' => 'dihedral_angle',
                    'units' => 'radians'
                };
            } else {
                $bond_parameters{$residue_name}{$dihedral_angle_name} = {
                    'values' =>
                        sample_angles( [ [ $dihedral_angle_start * $pi / 180.0,
                                           $dihedral_angle_end * $pi / 180.0 ] ],
                                       $dihedral_angle_count ),
                    'type' => 'dihedral_angle',
                    'units' => 'degrees'
                };
            }
        }
    }

    # Query overwrites on top.
    if( $query_strings ) {
        undef %bond_parameters;
    }

    my $residue_names_regexp = join '|', @{ $rotatable_residue_names };
    my $bond_parameter_regexp = '[A-Za-z0-9\-\.]+';
    my $float_regexp = '-?\d+(?:\.\d+)?';
    my $float_pos_regexp = '\d+(?:\.\d+)?';
    my $step_only_regexp =
        $legacy_grammar ? ${float_pos_regexp} : "\.\.(${float_pos_regexp})\.\.";

    for my $query_string ( split /;/, $query_strings ) {
        my $residue_names;
        my $bond_parameter_string;

        my @query_string_decomposed = split /:/, $query_string;
        if( scalar @query_string_decomposed == 2 ) {
            if( $query_string =~ m/^((?:${residue_names_regexp})(?:,(?:${residue_names_regexp}))*):(.+)$/i ) {
                $residue_names = [ split /,/, uc( $1 ) ];
                $bond_parameter_string = $2;
            } else {
                die "Syntax '$query_string' is incorrect\n";
            }
        } elsif( scalar @query_string_decomposed == 1 ) {
            $bond_parameter_string = $query_string;
        } else {
            die "Syntax '$query_string' is incorrect\n";
        }

        $residue_names //= [ "*" ];
        $bond_parameter_string //= "";

        for my $bond_parameter_query ( split /,/, $bond_parameter_string ) {
            my $bond_parameter_name;
            my $bond_parameter_start;
            my $bond_parameter_step;
            my $bond_parameter_end;

            if( $bond_parameter_query =~ m/^(${bond_parameter_regexp})=(${float_regexp})\.\.(${float_pos_regexp})\.\.(${float_regexp})$/ ) {
                ( $bond_parameter_name,
                  $bond_parameter_start,
                  $bond_parameter_step,
                  $bond_parameter_end ) = ( $1, $2, $3, $4 );
            } elsif( $bond_parameter_query =~ m/^(${bond_parameter_regexp})=(${float_regexp})\.\.(${float_regexp})$/ ) {
                ( $bond_parameter_name,
                  $bond_parameter_start,
                  $bond_parameter_end ) = ( $1, $2, $3 );
            } elsif( $bond_parameter_query =~ m/^(${bond_parameter_regexp})=${step_only_regexp}$/ ) {
                ( $bond_parameter_name, $bond_parameter_step ) = ( $1, $2 );
            } elsif( $bond_parameter_query =~ m/^(${bond_parameter_regexp})=(${float_regexp}|\!)$/ ) {
                if( $2 =~ m/^\!$/ ) {
                    # HACK: special value as all the steps have to be positive
                    # integer. It is set to -1 in order to make default values
                    # later.
                    ( $bond_parameter_name, $bond_parameter_step ) = ( $1, -1 );
                } else {
                    ( $bond_parameter_name,
                      $bond_parameter_start,
                      $bond_parameter_step,
                      $bond_parameter_end ) = ( $1, $2, 1.0, $2 + 1.0 );
                }
            } elsif( $bond_parameter_query =~ m/^(${float_regexp})\.\.(${float_pos_regexp})\.\.(${float_regexp})$/ ) {
                ( $bond_parameter_start,
                  $bond_parameter_step,
                  $bond_parameter_end ) = ( $1, $2, $3 );
            } elsif( $bond_parameter_query =~ m/^(${float_regexp})\.\.(${float_regexp})$/ ) {
                ( $bond_parameter_start, $bond_parameter_end ) = ( $1, $2 );
            } elsif( $bond_parameter_query =~ m/^${step_only_regexp}$/ ) {
                ( $bond_parameter_step ) = ( $1 );
            } elsif( $bond_parameter_query =~ m/^(${float_regexp}|\!)$/ ) {
                if( $1 =~ m/^\!$/ ) {
                    # HACK: special value as all the steps have to be positive
                    # integer. It is set to -1 in order to make default values
                    # later.
                    $bond_parameter_step = -1;
                } else {
                    ( $bond_parameter_start,
                      $bond_parameter_step,
                      $bond_parameter_end ) = ( $1, 1.0, $1 + 1.0 );
                }
            } else {
                die "Syntax '$bond_parameter_query' is incorrect\n";
            }

            $bond_parameter_start //= -180.0;
            $bond_parameter_step //= 36.0;
            $bond_parameter_end //= 180.0;
            $bond_parameter_name //= '*-*-*-*';

            my ( $bond_parameter_type ) =
                detect_bond_parameter_type( $bond_parameter_name );

            my $bond_parameter_count =
                int( ( $bond_parameter_end - $bond_parameter_start ) /
                     $bond_parameter_step );
            my $do_inclusive_end = 0;
            if( $bond_parameter_type eq 'bond_length' ||
                $bond_parameter_type eq 'bond_angle' ) {
                $bond_parameter_count++;
                $do_inclusive_end = 1;
            }

            my $bond_parameter_units = 'radians';
            if( $bond_parameter_type eq 'bond_length' ) {
                $bond_parameter_units = 'angstrom';
            } elsif( ! $in_radians ) {
                $bond_parameter_units = 'degrees';
                $bond_parameter_start = $bond_parameter_start * $pi / 180.0;
                $bond_parameter_end = $bond_parameter_end * $pi / 180.0;
            }

            for my $residue_name ( @{ $residue_names } ) {
                # HACK: not sure if the calculations above should be done if the
                # step value is -1.
                if( $bond_parameter_step < 0 ) {
                    $bond_parameters{$residue_name}{$bond_parameter_name} = {
                        'values' => [],
                        'type' => $bond_parameter_type,
                        'units' => $bond_parameter_units
                    };
                } else {
                    $bond_parameters{$residue_name}{$bond_parameter_name} = {
                        'values' => sample_bond_parameters(
                            [ [ $bond_parameter_start, $bond_parameter_end ] ],
                            $bond_parameter_count,
                            1,
                            $do_inclusive_end
                        ),
                        'type' => $bond_parameter_type,
                        'units' => $bond_parameter_units
                    };
                }
            }
        }
    }

    return \%bond_parameters;
}

sub retrieve_dihedral_angle_params
{
    my ( $dihedral_angle_restraints, $residue_name, $angle_name, $params ) = @_;

    my %params = ();
    my ( undef, undef, $any_angle_name ) =detect_bond_parameter_type( $angle_name );
    for my $param ( @{ $params } ) {
        if( exists $dihedral_angle_restraints->{$residue_name}{$angle_name}{$param} &&
            $dihedral_angle_restraints->{$residue_name}{$angle_name}{$param} ne '*' ) {
            $params{$param} =
                $dihedral_angle_restraints->{$residue_name}{$angle_name}{$param};
        } else {
            $params{$param} =
                $dihedral_angle_restraints->{'*'}{$any_angle_name}{$param};
        }
    }

    return map { $params{$_} } @{ $params };
}

1;
