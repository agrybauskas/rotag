%include LinearAlgebra.i

%module LinearAlgebra
%{
    #include "LinearAlgebra.h"
%}

%include "LinearAlgebra.h"

%perlcode %{
    use strict;
    use warnings;

    use Exporter qw( import );
    our @EXPORT_OK = qw( create_ref_frame );
%}
