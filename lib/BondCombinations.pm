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

    # Finds the longest bond parameter chain.
    my %parameter_keys =
        map { $_ => [ split /,/, $_ ]  } keys %{ $permuted_values };
    my @parameter_keys_sorted =
        sort { scalar( @{ $parameter_keys{$b} } ) <=>
               scalar( @{ $parameter_keys{$a} } ) ||
               $a cmp $b }
        keys %parameter_keys;

    # # Constructs bond parameter dependency tree.
    # my %parameter_tree = ();
    # for my $i ( 0..$#parameter_keys_sorted ) {
    #     my $parameter_key_1 = $parameter_keys_sorted[$i];
    #     for my $j ( $i+1..$#parameter_keys_sorted ) {
    #         my $parameter_key_2 = $parameter_keys_sorted[$j];
    #     }
    # }

    # Finds unique values for each parameter.
    my %visited_parameters = ();
    my %unique_parameter_values = ();
    while( @parameter_keys_sorted ) {
        my $parameter_key_sorted = shift @parameter_keys_sorted;
        my $parameter_names = $parameter_keys{$parameter_key_sorted};
        my $parameter_values = $permuted_values->{$parameter_key_sorted};

        next if ! grep { ! exists $visited_parameters{$_} } @{ $parameter_names };

        for my $i ( 0..$#{ $parameter_names } ) {
            my $parameter_name = $parameter_names->[$i];

            $visited_parameters{$parameter_name} = 1;

            my @unique_parameter_values =
                sort { $a <=> $b } uniq map { $_->[$i] } @{ $parameter_values };
            if( exists $unique_parameter_values{$parameter_name} ) {
                my %set1 =
                    map { $_ => 1 }
                    @{ $unique_parameter_values{$parameter_name} };
                my %set2 =
                    map { $_ => 1 }
                    @unique_parameter_values;
                $unique_parameter_values{$parameter_name} =
                    [ sort { $a <=> $b }
                      map { $_ }
                      grep { $set1{$_} }
                      keys %set2 ];
            } else {
                $unique_parameter_values{$parameter_name} =
                    clone \@unique_parameter_values;
            }
        }

        last if scalar( grep { exists $visited_parameters{$_} } @{ $names } ) ==
                scalar( @{ $names } );
    }

    # Generates permuted list from unique parameter values.
    my $combined_values = [];
    for my $parameter_name ( @{ $names } ) {
        if( @{ $combined_values } ) {
            $combined_values =
                permutation(
                    2, [],
                    [ $combined_values,
                      [ map { [ $_ ] }
                        @{ $unique_parameter_values{$parameter_name} } ] ], [] );
            $combined_values =
                [ map { [ @{ $_->[0] }, @{ $_->[1] } ] }
                     @{ $combined_values } ];
        } else {
            $combined_values =
                [ map { [ $_ ] }
                     @{ $unique_parameter_values{$parameter_name} } ];
        }
    }

    return $combined_values;
}

1;
