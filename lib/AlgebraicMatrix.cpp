#include "AlgebraicMatrix.h"

AlgebraicMatrix::AlgebraicMatrix(){};

/* -------------------------- Setters and Getters ---------------------------- */

void AlgebraicMatrix::set_symbols( std::vector<std::string> symbols )
{
  this->symbols = symbols;
  return;
}

void AlgebraicMatrix::set_is_evaluated( bool is_evaluated )
{
  this->is_evaluated = is_evaluated;
  return;
}

void AlgebraicMatrix::set_matrix( std::vector< std::vector<double > > matrix )
{
  this->matrix = matrix;
  return;
}

void AlgebraicMatrix::set_matrix_function( std::vector< std::vector<double> > ( *matrix_function )( double *args[] ) )
{
  this->matrix_function = matrix_function;
  return;
}
