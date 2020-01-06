#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::ParticleSwarm;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $particle_swarm =
    ParticleSwarm->new( { 'x' => { 'min' => -4.0, 'max' => 4.0 } }, 10 );

$particle_swarm->set_cost_function( \&{sub { my ( $x ) = @_; return $x**2 }} );
$particle_swarm->optimize( 10 );
