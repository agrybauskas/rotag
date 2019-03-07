#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "AlgebraicMatrix.h"
#include "LinearAlgebraCpp.h"

std::vector< std::vector<double> > matrix2_function( double args[] )
{
  double x = args[0];
  double y = args[1];
  double z = args[2];
  std::vector< std::vector<double> > tmp_matrix2 = { { x },
                                                     { y },
                                                     { z } };
  return tmp_matrix2;
}

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > matrix1 = { { 1, 0, 0 },
                                                 { 0, 2, 0 },
                                                 { 0, 0, 3 } };
  AlgebraicMatrix algebraic_matrix1;
  algebraic_matrix1.set_matrix( matrix1 );

  AlgebraicMatrix algebraic_matrix2;
  std::vector<std::string> symbols2 = { "x", "y", "z" };
  algebraic_matrix2.set_symbols( symbols2 );
  algebraic_matrix2.set_matrix_function( matrix2_function );

  std::map<std::string, double> symbol_values = { { "x", 6 },
                                                  { "y", 3 },
                                                  { "z", 2 } };

  std::vector<AlgebraicMatrix> matrices;
  matrices.push_back( algebraic_matrix1 );
  matrices.push_back( algebraic_matrix2 );

  mult_matrix_product( matrices, symbol_values );

  return 0;
}
