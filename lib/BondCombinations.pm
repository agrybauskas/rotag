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
    for my $name ( @{ $names } ) {
        for my $parameter_key ( @parameter_keys_sorted ) {
            if( $parameter_key =~ m/\Q$name\E/ ) {
                my @parameter_key_components = split /,/, $parameter_key;
                for my $i ( 0..$#parameter_key_components ) {
                    if( $name eq $parameter_key_components[$i] ) {
                        $parameter_name_pos{$name} = {
                            'pos' => $i,
                            'parameter_key' => $parameter_key,
                        }
                    }
                }
                last;
            }
        }
    }

    my $combined_values = [];
    for my $parameter_name ( @{ $names } ) {
        if( @{ $combined_values } ) {
            my $current_parameter_key =
                $parameter_name_pos{$parameter_name}{'parameter_key'};
            my $parameter_pos = $parameter_name_pos{$parameter_name}{'pos'};
            my @parameter_values =
                map { [ $_->[$parameter_pos] ] }
                   @{ $permuted_values->{$current_parameter_key} };
        } else {
            $combined_values =
                [ map { [ $_ ] }
                     @{ $permuted_values->{$parameter_name} } ];
        }
    }

    return $combined_values;
}

1;
