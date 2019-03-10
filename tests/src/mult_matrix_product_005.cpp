#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "AlgebraicMatrix.h"
#include "LinearAlgebra.h"

std::vector< std::vector<double> > matrix1_function( double args[] )
{
  double chi = args[0];
  std::vector< std::vector<double> > tmp_matrix1 = { { cos( chi ), -sin( chi ), 0 },
                                                     { sin( chi ),  cos( chi ), 0 },
                                                     { 0, 0, 1 } };
  return tmp_matrix1;
}

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
  std::vector<std::string> symbols1 = { "chi" };
  algebraic_matrix1.set_symbols( symbols1 );
  algebraic_matrix1.set_matrix_function( matrix1_function );

  AlgebraicMatrix algebraic_matrix2;
  std::vector<std::string> symbols2 = { "x", "y", "z" };
  algebraic_matrix2.set_symbols( symbols2 );
  algebraic_matrix2.set_matrix_function( matrix2_function );

  std::map<std::string, double> symbol_values = { { "x",   6 },
                                                  { "y",   3 },
                                                  { "z",   0 },
                                                  { "chi", 0 } };

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
