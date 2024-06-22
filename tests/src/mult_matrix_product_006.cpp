#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <map>

#include "AlgebraicMatrix.h"
#include "LinearAlgebra.h"

std::vector< std::vector<double> > matrix2_function( double args[] )
{
  double chi = args[0];
  std::vector< std::vector<double> > tmp_matrix2 = { { cos( chi ), -sin( chi ), 0 },
                                                     { sin( chi ),  cos( chi ), 0 },
                                                     { 0, 0, 1 } };
  return tmp_matrix2;
}

std::vector< std::vector<double> > matrix3_function( double args[] )
{
  double x = args[0];
  double y = args[1];
  double z = args[2];
  std::vector< std::vector<double> > tmp_matrix3 = { { x },
                                                     { y },
                                                     { z } };
  return tmp_matrix3;
}

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > matrix1 = { { 1, 2, 2 } };

  AlgebraicMatrix algebraic_matrix1;
  algebraic_matrix1.set_matrix( matrix1 );

  AlgebraicMatrix algebraic_matrix2;
  std::vector<std::string> symbols2 = { "chi" };
  algebraic_matrix2.set_symbols( symbols2 );
  algebraic_matrix2.set_matrix_function( matrix2_function );

  AlgebraicMatrix algebraic_matrix3;
  std::vector<std::string> symbols3 = { "x", "y", "z" };
  algebraic_matrix3.set_symbols( symbols3 );
  algebraic_matrix3.set_matrix_function( matrix3_function );

  std::map<std::string, double> symbol_values = { { "x",   6 },
                                                  { "y",   3 },
                                                  { "z",   2 },
                                                  { "chi", 0 } };

  std::vector<AlgebraicMatrix> matrices;
  matrices.push_back( algebraic_matrix1 );
  matrices.push_back( algebraic_matrix2 );
  matrices.push_back( algebraic_matrix3 );

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