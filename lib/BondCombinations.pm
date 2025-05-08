package BondCombinations;

use strict;
use warnings;

use Exporter qw( import );
BEGIN {
our @EXPORT_OK = qw( combine_permuted_values
                     parameter_key )
}

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

    my @combined_values;

    # Find the longest bond parameter chain.
    my %parameter_keys =
        map { $_ => [ split /,/, $_ ]  } keys %{ $permuted_values };
    my @parameter_keys_sorted =
        sort { scalar( @{ $parameter_keys{$b} } ) <=>
               scalar( @{ $parameter_keys{$a} } ) ||
               $a cmp $b }
        keys %parameter_keys;

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

            if( exists $unique_parameter_values{$parameter_name} ) {

            } else {
                $unique_parameter_values{$parameter_name} =
                    [ uniq map { $_->[$i] } @{ $parameter_values } ];
            }
        }

        last if scalar( grep { exists $visited_parameters{$_} } @{ $names } ) ==
                scalar( @{ $names } );
    }

    return \@combined_values;
}

1;
