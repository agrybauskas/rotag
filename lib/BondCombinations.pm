package BondCombinations;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( combine_permuted_values
                     parameter_key )
}

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

    my @combined_values;

    # Find the longest bond parameter chain.
    my %parameter_keys =
        map { $_ => [ split /,/, $_ ]  } keys %{ $permuted_values };
    my @parameter_keys_sorted =
        sort { scalar( @{ $parameter_keys{$b} } ) <=>
               scalar( @{ $parameter_keys{$a} } ) ||
               $a cmp $b }
        keys %parameter_keys;

    # Iterate through parameter name chains and retrieve values.
    my %visited_parameters = ();
    while( @parameter_keys_sorted ) {
        my $parameter_key_sorted = shift @parameter_keys_sorted;
        my $parameter_names = $parameter_keys{$parameter_key_sorted};
        for my $parameter_name ( @{ $parameter_names } ) {
            $visited_parameters{$parameter_name} = 1;
        }

        last if scalar( grep { exists $visited_parameters{$_} } @{ $names } ) ==
                scalar( @{ $names } );
    }
    return \@combined_values;
}

1;
