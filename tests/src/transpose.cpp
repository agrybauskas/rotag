#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include "LinearAlgebra.h"

int main( int argc, char* argv[] ) {
  std::vector< std::vector<double> > matrix;

  std::string line;
  std::ifstream fh( argv[1] );
  if( fh.is_open() ) {
    while( getline( fh, line ) ) {
      std::vector<std::string> line_split;
      std::vector<double> row;
      boost::split( line_split, line, boost::is_any_of(" ")  );
      for ( int i = 0; i < line_split.size(); i++ ) {
        if ( line_split[i] != "" ) {
          row.push_back( std::stod( line_split[i] ) );
        }
      }

      matrix.push_back( row );
    }
    fh.close();
  }

  std::vector< std::vector<double> > transposed_matrix = transpose( matrix );

  for( int i = 0; i < transposed_matrix.size(); i++ ) {
    for( int j = 0; j < transposed_matrix[i].size(); j++ ) {
      printf( "%.3f ", transposed_matrix[i][j] );
    }
    std::cout << std::endl;
  }

  return 0;
}
