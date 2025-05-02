#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use BondCombIter;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $bond_comb = BondCombIter->new( { 'chi1' => [ 1.0, 2.0 ],
                                     'chi2' => [ 3.0, 4.0, 5.0 ] } );

print Dumper $bond_comb->get_all_names();
print Dumper $bond_comb->get_all_values();
