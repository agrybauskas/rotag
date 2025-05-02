#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use BondCombinations;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $bond_comb = BondCombinations->new( {
    'order' => [ 'chi1' ],
    'collection' => {
        'chi1' => [ 1.0, 2.0 ]
    },
} );

print Dumper $bond_comb->get_all_names();
print Dumper $bond_comb->get_all_values();

$bond_comb->add( [ { 'chi2' => [ 3.0, 4.0, 5.0 ] } ] );

print Dumper $bond_comb->get_all_names();
print Dumper $bond_comb->get_all_values();
