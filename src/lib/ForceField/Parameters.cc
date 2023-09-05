#include "Parameters.h"

#include <cmath>
#include <iostream>
#include <memory>

extern "C" {
  #include "cif_compiler.h"
  #include "cif_options.h"
  #include "datablock.h"
}

#include "CIF.h"

Parameters::Parameters(char* parameter_file) {
  cif_option_t compiler_options = cif_option_default();
  cexception_t inner;
  CIF* parameters =
    new_cif_from_cif_file(parameter_file, compiler_options, &inner);
  datablock_dump(parameters->datablock_list);
}

Parameters::~Parameters() {};

double Parameters::epsilon() {
  double epsilon = 1.0;
  while((1.0 + 0.5 * epsilon) != 1.0) {
      epsilon = 0.5 * epsilon;
  }
  return epsilon;
}

double Parameters::pi() {
  return 4 * std::atan2(1, 1);
}
