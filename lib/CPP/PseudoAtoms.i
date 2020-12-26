%module "PseudoAtoms"

%inline %{
    extern void generate_library(SV* hash_ref);
%}
