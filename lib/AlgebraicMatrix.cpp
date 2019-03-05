#include "AlgebraicMatrix.h"

AlgebraicMatrix::AlgebraicMatrix( args arguments )
{
  this->is_evaluated = arguments.is_evaluated;
  this->symbols = arguments.symbols;
  this->matrix = arguments.matrix;
  return;
}

void AlgebraicMatrix::set_symbols( std::vector<std::string> symbols )
{
  this->symbols = symbols;
  return;
}
