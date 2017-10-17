%module LinearAlgebra
%{
  #include "LinearAlgebra.h"
%}

%include "std_vector.i"
%template(VectorDouble) std::vector<double>;
%template(VectorVectorDouble) std::vector< std::vector<double> >;

%include "LinearAlgebra.h"
