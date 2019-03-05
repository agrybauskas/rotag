#include "AlgebraicMatrix.h"

AlgebraicMatrix::AlgebraicMatrix( bool is_evaluated,
                                  std::vector<std::string> symbols,
                                  std::vector< std::vector<double> > matrix,
                                  std::vector< std::vector<double> > ( *matrix_func )( double *args[] ) )
{
  this->is_evaluated = is_evaluated;
  this->symbols = symbols;
  this->matrix = matrix;
  return;
}

void AlgebraicMatrix::set_symbols( std::vector<std::string> symbols )
{
  this->symbols = symbols;
  return;
}
