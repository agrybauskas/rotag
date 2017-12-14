%include GinacAlgebra.i

%module GinacAlgebra
%{
  #include "GinacAlgebra.h"
%}

%include "std_vector.i"
%include "std_string.i"

%template(VecString) std::vector<string>;

%include "GinacAlgebra.h"

%perlcode %{
  use strict;
  use warnings;

  use Exporter qw( import );
  our @EXPORT_OK = qw( matrix_product_ginac );
%}
