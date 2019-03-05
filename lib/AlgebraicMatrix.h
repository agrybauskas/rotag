#ifndef _ALGEBRAICMATRIX_H_
#define _ALGEBRAICMATRIX_H_

#include <string>
#include <vector>

class AlgebraicMatrix {
 private:
  bool is_evaluated;
  std::vector<std::string> symbols;
  std::vector< std::vector<double> > matrix;
  std::vector< std::vector<double> > matrix_function;

 public:
  AlgebraicMatrix();
  void set_symbols( std::vector<std::string> symbols );
};

#endif
