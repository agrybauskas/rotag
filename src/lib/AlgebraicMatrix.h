#ifndef SRC_LIB_ALGEBRAICMATRIX_H_
#define SRC_LIB_ALGEBRAICMATRIX_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

class AlgebraicMatrix {
 private:
  bool is_evaluated = false;
  std::vector<std::string> symbols;
  std::vector< std::vector<double>> matrix;
  std::vector< std::vector<double>>(*matrix_function)(double args[]);

 public:
  AlgebraicMatrix();

  void set_symbols(std::vector<std::string> symbols);
  void set_matrix(std::vector<std::vector<double>> matrix);
  void set_matrix_function(std::vector<std::vector<double>>
                           (*matrix_function)(double args[]));

  void evaluate(std::map<std::string, double> symbol_values);
  std::vector<std::string> get_symbols();
  bool get_is_evaluated();
  std::vector< std::vector<double > > get_matrix();
};

#endif  // SRC_LIB_ALGEBRAICMATRIX_H_