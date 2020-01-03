#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::Parameter;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $param = Parameter->new( { 'key' => 'epsilon',
                              'min_range' => 3.0,
                              'max_range' => 5.0 } );

print Dumper $param->key;
print Dumper $param->min_range;
print Dumper $param->max_range;
