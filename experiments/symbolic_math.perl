use Math::Algebra::Symbols;
use Math::Complex;
use Data::Dumper;
  
use strict;
use warnings;

print( "-" x 80, "\n" );

#
# Trying cos(x)^2 + sin(x)^2 = 1 type expression.
#

my ( $chi ) = symbols( qw( chi ) );

my $expression = $chi**2;

$chi = 2;

print( eval( $expression ), "\n" );

print( "-" x 80, "\n" );
