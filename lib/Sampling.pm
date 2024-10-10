package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_bond_parameters
                     sample_bond_parameters_qs_parsing );

use POSIX;

use BondParameters qw( detect_bond_parameter_type );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Sampling ---------------------------------- #

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

    my $rotatable_residue_names =
        $parameters->{'_[local]_rotatable_residue_names'};
    my $bond_parameter_restraints =
        $parameters->{'_[local]_bond_parameter_restraints'};

    $query_strings =~ s/\s//g;

    # Dihedral angle calculations are first preditermined and then updated.
    my %bond_parameters;

    # Default bond parameters.
    for my $residue_name ( sort keys %{ $bond_parameter_restraints } ) {
        for my $bond_parameter_name (
            sort keys %{ $bond_parameter_restraints->{$residue_name} } ) {
            my ( $bond_parameter_start,
                 $bond_parameter_step,
                 $bond_parameter_end ) =
                map { $bond_parameter_restraints->{$residue_name}
                                                  {$bond_parameter_name}{$_} }
                    ( 'range_from', 'step', 'range_to' );
            my ( $bond_parameter_type ) =
                detect_bond_parameter_type( $bond_parameter_name );
            my $bond_parameter_values = determine_bond_parameter_values(
                $parameters,
                $bond_parameter_name,
                $bond_parameter_start,
                $bond_parameter_step,
                $bond_parameter_end,
                { 'in_radians' => $in_radians }
            );
            $bond_parameters{$residue_name}{$bond_parameter_name} = {
                'from' => $bond_parameter_start,
                'step' => $bond_parameter_step,
                'to' => $bond_parameter_end,
                'values' => $bond_parameter_values,
                'type' => $bond_parameter_type,
                'units' => (
                    $bond_parameter_type eq 'bond_length' ?
                    'angstroms' :
                    ( $in_radians ? 'radians' : 'degrees' )
                )
            };
        }
    }

    # Overwrites the parameters if it was declared.
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

            $bond_parameter_name //= '*-*-*-*';

            my $bond_parameter_values = determine_bond_parameter_values(
                $parameters,
                $bond_parameter_name,
                $bond_parameter_start,
                $bond_parameter_step,
                $bond_parameter_end,
                { 'in_radians' => $in_radians }
            );

            for my $residue_name ( @{ $residue_names } ) {
                $bond_parameters{$residue_name}{$bond_parameter_name} = {
                    'from' => $bond_parameter_start,
                    'step' => $bond_parameter_step,
                    'to' => $bond_parameter_end,
                    'values' => (
                        $bond_parameter_step < 0 ? [] : $bond_parameter_values
                    ),
                    'type' => detect_bond_parameter_type( $bond_parameter_name ),
                    'units' => ( $in_radians ? 'radians' : 'degrees' )
                };
            }
        }
    }

    resolve_bond_parameters( \%bond_parameters );

    return \%bond_parameters;
}

sub determine_bond_parameter_values
{
    my ( $parameters, $bond_parameter_name, $bond_parameter_start,
         $bond_parameter_step, $bond_parameter_end, $options ) = @_;

    return if $bond_parameter_start eq '*' ||
              $bond_parameter_step eq '*' ||
              $bond_parameter_end eq '*';

    my ( $in_radians ) = ( $options->{'in_radians'} );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};

    my ( $bond_parameter_type ) =
        detect_bond_parameter_type( $bond_parameter_name );
    my $bond_parameter_count =
        int( ( $bond_parameter_end - $bond_parameter_start ) /
             $bond_parameter_step );

    my $values;
    if( $bond_parameter_type eq 'dihedral_angle' && $in_radians ) {
        $values = sample_bond_parameters(
            [ [ $bond_parameter_start,
                $bond_parameter_end ] ],
            $bond_parameter_count, 1, 0
        );
    } elsif( $bond_parameter_type eq 'dihedral_angle' ) {
        $values = sample_bond_parameters(
            [ [ $bond_parameter_start * $pi / 180.0,
                $bond_parameter_end * $pi / 180.0 ] ],
            $bond_parameter_count, 1, 0
        );
    } else {
        $values = sample_bond_parameters(
            [ [ $bond_parameter_start,
                $bond_parameter_end ] ],
            $bond_parameter_count, 1, 1
        );
    }
    return $values;
}

sub resolve_bond_parameters
{
    my ( $bond_parameters ) = @_;
    my %resolved_bond_parameters;
    for my $residue_name ( keys %{ $bond_parameters } ) {
        for my $bond_parameter ( keys %{ $bond_parameters->{$residue_name} } ) {
        }
    }

    return;
}

1;
