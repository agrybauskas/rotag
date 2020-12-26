%module "PseudoAtoms"

%inline %{
    extern void generate_library(SV* options_ptr);
%}
