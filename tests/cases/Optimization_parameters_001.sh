#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::Parameters;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $optimization_param = Parameters->new( { 'key' => 'epsilon',
                                            'min_range' => 3.0,
                                            'max_range' => 5.0,
                                            'value' => 5.0 } );

$optimization_param->speed( 4 );

print Dumper $optimization_param->key;
print Dumper $optimization_param->min_range;
print Dumper $optimization_param->max_range;
print Dumper $optimization_param->speed;
print Dumper $optimization_param->value;
