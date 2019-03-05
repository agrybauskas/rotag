#include "AlgebraicMatrix.h"

AlgebraicMatrix::AlgebraicMatrix(){};

void AlgebraicMatrix::set_symbols( std::vector<std::string> symbols )
{
  this->symbols = symbols;
  return;
}
