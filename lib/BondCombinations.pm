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
    my @parameter_keys_sorted_rev = reverse @parameter_keys_sorted;

    # Constructs bond parameter one way dependency tree.
    my %parameter_tree = ();
    for my $i ( 0..$#parameter_keys_sorted_rev ) {
        my $parameter_key_1 = $parameter_keys_sorted_rev[$i];
        for my $j ( $i+1..$#parameter_keys_sorted_rev ) {
            my $parameter_key_2 = $parameter_keys_sorted_rev[$j];
            if( $parameter_key_2 =~ m/^\Q$parameter_key_1\E/ ) {
                $parameter_tree{$parameter_key_2} = $parameter_key_1;
            }
        }
    }

    # Identifies the longest parameter key for given parameter name.
    my %parameter_longest = ();
    for my $parameter_name ( @{ $names } ) {
        for my $parameter_key ( @parameter_keys_sorted_rev ) {
            if( $parameter_key =~ m/\Q$parameter_name\E$/ ) {
                $parameter_longest{$parameter_name} = $parameter_key;
                last;
            }
        }
    }

    # Generates permuted list from unique parameter values.
    my $combined_values = [];
    for my $parameter_name ( @{ $names } ) {
        if( @{ $combined_values } ) {
        } else {
            $combined_values =
                [ map { [ $_ ] }
                     @{ $permuted_values->{$parameter_name} } ];
        }
    }

    return $combined_values;
}

1;
