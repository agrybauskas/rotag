#include "Parameters.h"

extern "C" {
  #include "cif_compiler.h"
  #include "cif_options.h"
}

Parameters::Parameters(char* parameter_file) {
    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    CIF* parameters =
        new_cif_from_cif_file(parameter_file, compiler_options, &inner);
}
