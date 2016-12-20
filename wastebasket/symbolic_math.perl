use Math::Algebra::Symbols;
  
my $theta = symbols( qw( theta ) );
my $expression = 2 * $theta + $theta;
my $theta = 5;

print eval( $expression );
