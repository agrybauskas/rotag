%include LinearAlgebraCpp.i

%module(package="LinearAlgebraCpp") "LinearAlgebraCpp"
%{
    #include "LinearAlgebraCpp.h"
%}

%include "LinearAlgebraCpp.h"

%perlcode %{
    use strict;
    use warnings;

    use Exporter qw( import );
    our @EXPORT_OK = qw( create_ref_frame );
%}
