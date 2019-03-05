#ifndef _ALGEBRAICMATRIX_H_
#define _ALGEBRAICMATRIX_H_

#include <string>
#include <vector>

class AlgebraicMatrix {
 private:
  bool is_evaluated;
  std::vector<std::string> symbols;
  std::vector< std::vector<double> > matrix;

 public:
  AlgebraicMatrix( bool is_evaluated,
                   std::vector<std::string> symbols,
                   std::vector<std::string> *matrix_func() );

  void set_symbols( std::vector<std::string> symbols );
};

#endif
