#!/usr/bin/perl

use Math::Algebra::Symbols;
use Math::Complex;
use Data::Dumper;

use strict;
use warnings;

#
# Trying cos(x)^2 + sin(x)^2 = 1 type expression.
#

my $PI = 4 * atan2( 1, 1 );
my ( $chi ) = symbols( qw( chi ) );
my $expression = cos( $chi ) + cos( $chi );

print( $expression, "\n");
print( $expression->{"s"}, "\n" );
