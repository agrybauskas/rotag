package BondCombinations;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( combine_permuted_values
                     parameter_key )
}

use Clone qw( clone );
use Combinatorics qw( permutation );
use List::Util qw( uniq );

sub parameter_key
{
    my ( $names ) = @_;
    return join ',', @{ $names };
}

sub combine_permuted_values
{
    my ( $permuted_values, $names ) = @_;
    my $parameter_key = parameter_key( $names );

    return $permuted_values->{$parameter_key}
        if exists $permuted_values->{$parameter_key};

    my %parameter_keys =
        map { $_ => [ split /,/, $_ ]  } keys %{ $permuted_values };
    my @parameter_keys_sorted =
        sort { scalar( @{ $parameter_keys{$b} } ) <=>
               scalar( @{ $parameter_keys{$a} } ) ||
               $a cmp $b }
        keys %parameter_keys;

    my %parameter_name_pos = ();
    my %parameter_name_groups = ();
    for my $name ( @{ $names } ) {
        for my $parameter_key ( @parameter_keys_sorted ) {
            if( $parameter_key =~ m/\Q$name\E/ ) {
                my @parameter_key_components = split /,/, $parameter_key;
                for my $i ( 0..$#parameter_key_components ) {
                    if( $name eq $parameter_key_components[$i] ) {
                        $parameter_name_pos{$name} = {
                            'pos' => $i,
                            'parameter_key' => $parameter_key,
                        };
                        push @{ $parameter_name_groups{$parameter_key} },
                            $name;
                    }
                }
                last;
            }
        }
    }

    my %visited_parameter_keys = ();
    my $combined_values = [];
    for my $parameter_name ( @{ $names } ) {
        my $parameter_key =
            $parameter_name_pos{$parameter_name}{'parameter_key'};

        next if $visited_parameter_keys{$parameter_key};
        $visited_parameter_keys{$parameter_key} = 1;

        my @parameter_name_comb = @{ $parameter_name_groups{$parameter_key} };
        my @parameter_value_comb = ();
        for my $current_parameter_name ( @parameter_name_comb ) {
            my $parameter_pos =
                $parameter_name_pos{$current_parameter_name}{'pos'};
            push @parameter_value_comb,
                [ map { [ $_->[$parameter_pos] ] }
                     @{ $permuted_values->{$parameter_key} } ];
        }

        my @parameter_value_comb_flattened = ();
        for my $param_idx ( 0..$#parameter_value_comb ) {
            for my $i ( 0..$#{ $parameter_value_comb[$param_idx] } ) {
                if( $param_idx > 0 ) {
                    push @{ $parameter_value_comb_flattened[$i] },
                        @{ $parameter_value_comb[$param_idx][$i] };
                } else {
                    push @parameter_value_comb_flattened,
                        $parameter_value_comb[$param_idx][$i];
                }
            }
        }

        if( @{ $combined_values } ) {
            $combined_values =
                permutation(
                    2, [],
                    [ $combined_values,
                      \@parameter_value_comb_flattened ], [] );
            $combined_values =
                [ map { [ @{ $_->[0] }, @{ $_->[1] } ] }
                     @{ $combined_values } ];
        } else {
            $combined_values = [ @parameter_value_comb_flattened ];
        }
    }

    # Value positions have to be reordered as bond parameter keys might vary
    # positionally.
    my @combined_values_reordered = ();
    for my $name ( @{ $names } ) {

    }

    return \@combined_values_reordered;
}

1;
