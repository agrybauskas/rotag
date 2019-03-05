#include "AlgebraicMatrix.h"

#include <iostream>
#include <string>
#include <vector>

std::vector< std::vector<double> > matrix2_function( double *args[] )
{
    std::vector< std::vector<double> > tmp_matrix2 = { { 1, 0, 0 },
                                                       { 0, 2, 0 },
                                                       { 0, 0, 3 } };
    return tmp_matrix2;
}

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > matrix1 = { { 1, 0, 0 },
                                                 { 0, 2, 0 },
                                                 { 0, 0, 3 } };
  AlgebraicMatrix algebraic_matrix1;
  algebraic_matrix1.set_matrix( matrix1 );

  std::cout << algebraic_matrix1.get_is_evaluated() << std::endl;

  AlgebraicMatrix algebraic_matrix2;
  std::vector<std::string> symbols2 = { "x", "y", "z" };
  algebraic_matrix2.set_symbols( symbols2 );
  algebraic_matrix2.set_matrix_function( matrix2_function );

  std::cout << algebraic_matrix2.get_is_evaluated() << std::endl;

  return 0;
}
