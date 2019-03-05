#include "AlgebraicMatrix.h"

#include <functional>
#include <string>
#include <vector>

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > tmp_matrix1 = { { 1, 0, 0 },
                                                     { 0, 2, 0 },
                                                     { 0, 0, 3 } };

  AlgebraicMatrix matrix1( true, std::vector<std::string> (), tmp_matrix1 );

  std::vector<std::string> symbols2 = { "x", "y", "z" };
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
