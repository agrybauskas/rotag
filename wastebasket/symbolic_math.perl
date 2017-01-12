use Math::Algebra::Symbols;
use Data::Dumper;
  
use strict;
use warnings;

print( "-" x 80, "\n" );

#
# Trying cos(x)^2 + sin(x)^2 = 1 type expression.
#

my ( $chi ) = symbols( qw( chi ) );

my $expression = ( cos($chi)**2 ) + ( sin($chi)**2 );

$chi = 4;

print( eval( '$expression' ), "\n" );

print( "-" x 80, "\n" );
