#include "AlgebraicMatrix.h"

/* --------------------- Contructors and/or destructors ---------------------- */

AlgebraicMatrix::AlgebraicMatrix(){};

/* -------------------------------- Methods ---------------------------------- */

void AlgebraicMatrix::evaluate( std::map<std::string, double> symbol_values )
{
  std::vector<std::string> symbols = this->symbols;

  /* Checks if all symbol values are present. */
  try {
    for( std::string symbol : symbols ) {
      if( symbol_values.count( symbol ) == 0 ) {
        throw "'" + symbol + "'" + " value is not passed";
      }
    }
  } catch( std::string error_message ) {
    std::cout << error_message << std::endl;
    exit( EXIT_FAILURE );
  }

  /* Prepares array for passing to matrix function. */
  int symbols_size = symbols.size();
  double args[symbols_size];
  for( int i = 0; i < symbols_size; i++ ) {
    args[i] = symbol_values[symbols[i]];
  }

  this->matrix = this->matrix_function( args );
  this->is_evaluated = true;

  return;
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
    std::vector< std::vector<double> > ( *matrix_function )( double args[] ) )
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
