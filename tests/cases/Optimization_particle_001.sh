#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::Particle;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $particle = Particle->new( { 'epsilon' => { 'min' => 3.0, 'max' => 5.0 } } );

$particle->position( [ 4, 5, 6 ] );
$particle->speed( [ -4, -5, 6 ] );

print Dumper $particle->position;
print Dumper $particle->speed;
