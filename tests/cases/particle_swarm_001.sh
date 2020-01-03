#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Optimization qw( particle_swarm );

# $Data::Dumper::Sortkeys = 1;
# $Data::Dumper::Indent = 1;

# my $particles =
#     Optimization->new( { 'x'   => { 'min_range' => 1.0,
#                                     'max_range' => 2.0 } }, 10 );
# $particles->set_cost_function( \&{ sub { my ( $x ) = @_; return $x^2 } } );

# particle_swarm( $particles );
