package Sampling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sample_angles
                     sample_bond_parameters
                     sample_angles_qs_parsing
                     sample_bond_parameters_qs_parsing );

use POSIX;

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

            my $angle_count = int( ( $angle_end - $angle_start ) / $angle_step );

            if( $in_radians ) {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( [ [ $angle_start, $angle_end ] ],
                                   $angle_count );
            } else {
                $angles{$residue_name}{$angle_name} =
                    sample_angles( [ [ $angle_start * $pi / 180.0,
                                       $angle_end * $pi / 180.0 ] ],
                                   $angle_count );
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

            my $angle_count = int( ( $angle_end - $angle_start ) / $angle_step );

            for my $residue_name ( @{ $residue_names } ) {
                if( $in_radians ) {
                    $angles{$residue_name}{$angle_name} =
                        sample_angles( [ [ $angle_start, $angle_end ] ],
                                       $angle_count );
                } else {
                    $angles{$residue_name}{$angle_name} =
                        sample_angles( [ [ $angle_start * $pi / 180.0,
                                           $angle_end * $pi / 180.0 ] ],
                                       $angle_count );
                }
            }
        }
    }

    return \%angles;
}

sub sample_bond_parameters_qs_parsing
{
    my ( $parameters, $query_strings, $options ) = @_;

    my ( $in_radians ) = ( $options->{'in_radians'} );

    my $pi = $parameters->{'_[local]_constants'}{'pi'};
    my $dihedral_angle_restraints =
        $parameters->{'_[local]_dihedral_angle_restraints'};
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
            if( $residue_name eq '.' ) {
                $residue_name = '*'
            }
            if( $dihedral_angle_name eq '.' ) {
                $dihedral_angle_name = '*'
            }
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
                $bond_parameters{$residue_name}
                                {'rotatable_bonds'}
                                {$dihedral_angle_name} =
                    sample_angles( [ [ $dihedral_angle_start,
                                       $dihedral_angle_end ] ],
                                   $dihedral_angle_count );
            } else {
                $bond_parameters{$residue_name}
                                {'rotatable_bonds'}
                                {$dihedral_angle_name} =
                    sample_angles( [ [ $dihedral_angle_start * $pi / 180.0,
                                       $dihedral_angle_end * $pi / 180.0 ] ],
                                   $dihedral_angle_count );
            }
        }
    }

    # Query overwrites on top.
    if( $query_strings ) {
        undef %bond_parameters;
    }

    for my $query_string ( split /;/, $query_strings ) {
        my $residue_names;
        my $bond_parameter_string;

        my $residue_names_regexp = join '|', @{ $rotatable_residue_names };
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

        for my $bond_parameter ( split /,/, $bond_parameter_string ) {
            my $bond_parameter_name;
            my $bond_parameter_start;
            my $bond_parameter_step;
            my $bond_parameter_end;

            if( $bond_parameter =~ m/^([A-Za-z0-9\-]+)=(-?\d+(?:\.\d+)?)\.\.(\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_name,
                  $bond_parameter_start,
                  $bond_parameter_step,
                  $bond_parameter_end ) = ( $1, $2, $3, $4 );
            } elsif( $bond_parameter =~ m/^([A-Za-z0-9\-]+)=(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_name,
                  $bond_parameter_start,
                  $bond_parameter_end ) = ( $1, $2, $3 );
            } elsif( $bond_parameter =~ m/^([A-Za-z0-9\-]+)=(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_name, $bond_parameter_step ) = ( $1, $2 );
            } elsif( $bond_parameter =~ m/^(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_step ) = ( $1 );
            } elsif( $bond_parameter =~ m/^(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_start,
                  $bond_parameter_step,
                  $bond_parameter_end ) = ( $1, $2, $3 );
            } elsif( $bond_parameter =~ m/^(-?\d+(?:\.\d+)?)\.\.(-?\d+(?:\.\d+)?)$/ ) {
                ( $bond_parameter_start, $bond_parameter_end ) = ( $1, $2 );
            } else {
                die "Syntax '$bond_parameter' is incorrect\n";
            }

            $bond_parameter_name //= '*';
            $bond_parameter_start //= - 180.0;
            $bond_parameter_step //= 36.0;
            $bond_parameter_end //= 180.0;

            # Determine bond parameter type.
            # HACK: later, probably should also check force field
            # parameter file.
            my $bond_parameter_type;
            if( scalar( split( /-/, $bond_parameter_name ) ) == 3 ) {
                $bond_parameter_type = 'rotatable_bonds';
            } elsif( scalar( split( /-/, $bond_parameter_name ) ) == 2 ) {
                $bond_parameter_type = 'rotatable_bonds';
            } else {
                $bond_parameter_type = 'rotatable_bonds';
            }

            my $bond_parameter_count =
                int( ( $bond_parameter_end - $bond_parameter_start ) /
                     $bond_parameter_step );

            for my $residue_name ( @{ $residue_names } ) {
                if( $in_radians ) {
                    $bond_parameters{$residue_name}
                                    {$bond_parameter_type}
                                    {$bond_parameter_name} =
                        sample_angles( [ [ $bond_parameter_start,
                                           $bond_parameter_end ] ],
                                       $bond_parameter_count );
                } else {
                    $bond_parameters{$residue_name}
                                    {$bond_parameter_type}
                                    {$bond_parameter_name} =
                        sample_angles( [ [ $bond_parameter_start * $pi / 180.0,
                                           $bond_parameter_end * $pi / 180.0 ] ],
                                       $bond_parameter_count );
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
