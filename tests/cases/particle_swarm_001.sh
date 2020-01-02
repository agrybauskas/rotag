#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization qw( particle_swarm );

my $particles =
    Optimization->new( { 'phase'   => { 'min_range' => 1.0,
                                        'max_range' => 2.0 },
                         'epsilon' => { 'min_range' => 1.0,
                                        'max_range' => 5.0 } } );
particle_swarm( $particles );
