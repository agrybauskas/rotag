%module LinearAlgebra
%{
  #include "LinearAlgebra.h"
%}

%include "std_vector.i"

namespace std {
  %template(vector_double) vector<double>;
  %template(vector_vector_double) vector< vector<double> >;
}

%include "LinearAlgebra.h"
