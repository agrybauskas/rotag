#include "Parameters.h"

#include "cif_compiler.h"
#include "cif_options.h"

Parameters::Parameters(std::string parameter_file) {
    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    CIF* parameters =
        new_cif_from_cif_file((char*) &parameter_file, compiler_options, &inner);
}
