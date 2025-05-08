package BondCombinations;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( combine_permuted_values
                     parameter_key )
}

use List::Util qw( any );

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

    return \@combined_values;
}

1;
