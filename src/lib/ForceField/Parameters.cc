#include "Parameters.h"

#include "cif_compiler.h"

Parameters::Parameters(std::string parameter_file) {
    CIF* parameters = new_cif_from_cif_file((char*) &parameter_file);
}
