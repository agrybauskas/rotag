%include LinearAlgebra.i

%module LinearAlgebra
%{
  #include "LinearAlgebra.h"
%}

%include "std_vector.i"

%template(VecInt) std::vector<int>;
%template(VecVecInt) std::vector<std::vector<int> >;

%include "LinearAlgebra.h"

%perlcode %{
  use strict;
  use warnings;

  use Exporter qw( import );
  our @EXPORT_OK = qw( matrix_product );
%}
