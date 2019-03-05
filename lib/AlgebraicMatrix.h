#ifndef _ALGEBRAICMATRIX_H_
#define _ALGEBRAICMATRIX_H_

#include <string>
#include <vector>

struct args {
  bool is_evaluated;
  std::vector<std::string> symbols;
  std::vector< std::vector<double> > matrix;
};

class AlgebraicMatrix {
 private:
  bool is_evaluated;
  std::vector<std::string> symbols;
  std::vector< std::vector<double> > matrix;

 public:
  AlgebraicMatrix( args arguments );

  void set_symbols( std::vector<std::string> symbols );
};

#endif
