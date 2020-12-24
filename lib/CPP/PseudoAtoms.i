%include PseudoAtoms.i

%module(package="PseudoAtoms") "PseudoAtoms"
%{
    #include "PseudoAtoms.h"
%}

%include "PseudoAtoms.h"

%perlcode %{
    use strict;
    use warnings;

    use Exporter qw( import );
    our @EXPORT_OK = qw( generate_library );
%}
