#include "AlgebraicMatrix.h"

#include <functional>
#include <string>
#include <vector>

std::vector< std::vector<double> > matrix2_func( double *args[] )
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
    // AlgebraicMatrix matrix2( true, std::vector<std::string> (), matrix2_func );

  // std::vector<std::string> symbols2 = { "x", "y", "z" };
  // AlgebraicMatrix matrix2(
  //     false,
  //     symbols2,
  //     std::function< std::vector< std::vector<double> >() >
  //     {
  //       std::vector< std::vector<double> > tmp_matrix1 = { { 1, 0, 0 },
  //                                                          { 0, 2, 0 },
  //                                                          { 0, 0, 3 } };
  //       return tmp_matrix1;
  //     }
  // );

  return 0;
}
