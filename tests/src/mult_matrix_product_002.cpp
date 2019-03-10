#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "AlgebraicMatrix.h"
#include "LinearAlgebra.h"

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > matrix1 = { { 1, 0, 0 },
                                                 { 0, 2, 0 },
                                                 { 0, 0, 3 } };
  std::vector< std::vector<double> > matrix2 = { { 3 }, { 2 }, { 1 } };

  AlgebraicMatrix algebraic_matrix1;
  algebraic_matrix1.set_matrix( matrix1 );

  AlgebraicMatrix algebraic_matrix2;
  algebraic_matrix2.set_matrix( matrix2 );

  std::map<std::string, double> symbol_values;

  std::vector<AlgebraicMatrix> matrices;
  matrices.push_back( algebraic_matrix1 );
  matrices.push_back( algebraic_matrix2 );

  std::vector<std::vector<double>> matrix_product_result =
    mult_matrix_product( matrices, symbol_values )[0].get_matrix();

  for( int i = 0; i < matrix_product_result.size(); i++ ) {
    for( int j = 0; j < matrix_product_result[i].size(); j++ ) {
      printf( "%.3f ", matrix_product_result[i][j] );
    }
    std::cout << std::endl;
  }

  return 0;
}
