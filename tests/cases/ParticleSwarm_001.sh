#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization::ParticleSwarm;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $particle_swarm = ParticleSwarm->new();
