#include "AlgebraicMatrix.h"

/* --------------------- Contructors and/or destructors ---------------------- */

AlgebraicMatrix::AlgebraicMatrix(){};

/* -------------------------------- Methods ---------------------------------- */

void AlgebraicMatrix::evaluate( std::map<std::string, double> )
{

}

/* -------------------------- Setters and Getters ---------------------------- */

void AlgebraicMatrix::set_symbols( std::vector<std::string> symbols )
{
  this->symbols = symbols;
  return;
}

void AlgebraicMatrix::set_matrix( std::vector< std::vector<double > > matrix )
{
  this->matrix = matrix;
  this->is_evaluated = true;
  return;
}

void AlgebraicMatrix::set_matrix_function(
    std::vector< std::vector<double> > ( *matrix_function )( double *args[] ) )
{
  this->matrix_function = matrix_function;
  return;
}

std::vector<std::string> AlgebraicMatrix::get_symbols() {
  return this->symbols;
}

bool AlgebraicMatrix::get_is_evaluated()
{
  return this->is_evaluated;
}

std::vector< std::vector<double > > AlgebraicMatrix::get_matrix()
{
  return this->matrix;
}
