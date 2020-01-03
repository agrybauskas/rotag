#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::Particle;
use Optimization qw( particle_swarm );

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $x = Particle->new( { 'x' => { 'min_range' => 5.0, 'max_range' => 10.0 } } );

particle_swarm( [ $x ] );
