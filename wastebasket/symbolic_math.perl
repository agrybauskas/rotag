use Math::Algebra::Symbols;
use Data::Dumper;
  
my ( $theta, $chi ) = symbols( qw( theta chi ) );

#
# Simple variable substitution.
#

my $expression = 2 * $theta + $theta;

print( $expression, "\n" );

print( "-" x 80, "\n" );

#
# Variable behaviour inside arrays.
#

my @matrices = ( ( $chi ), ( $chi * 2 ) );

foreach( @matrices ) {
    print( $_, "\n" );
}

print( "-" x 80, "\n" );

#
# Combining multiple arrays into one expression.
#

my @matrices = ( ( ( $chi + $chi, $theta, 1 ),
		   ( $theta, $chi, 0 ) ),

		 ( ( $theta, $chi, 3 ),
		   ( 0,      1,    0 ) ) );

my $sum_of_all = 0;

foreach( @matrices ) {
    foreach( $_ ) {
	$sum_of_all += $_;
    }
}

print( $sum_of_all, "\n" );

print( "-" x 80, "\n" );

#
# Assigning symbols by argvar.
#

my @symbols = qw( $theta $lambda );

foreach( @symbols ) {
    eval( "$_ = 25" );
}

print( $theta + $lambda, "\n" );

print( "-" x 80, "\n" );

print Dumper $chi*$chi;

print( "-" x 80, "\n" );
