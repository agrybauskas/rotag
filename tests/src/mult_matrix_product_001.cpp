#include "AlgebraicMatrix.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

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

  algebraic_matrix2.evaluate( symbol_values );

  for( int i = 0; i < algebraic_matrix2.get_matrix().size(); i++ ) {
    for( int j = 0; j < algebraic_matrix2.get_matrix()[i].size(); j++ ) {
      printf( "%.3f ", algebraic_matrix2.get_matrix()[i][j] );
    }
    std::cout << std::endl;
  }

  return 0;
}
