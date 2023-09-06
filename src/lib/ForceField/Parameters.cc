#include "Parameters.h"

#include <cmath>
#include <iostream>
#include <memory>

#include <boost/filesystem.hpp>

extern "C" {
  #include "cif_compiler.h"
}

Parameters::Parameters() {
  // boost::filesystem::path parameter_file_path{parameter_file};
  // std::cout << boost::filesystem::current_path() << std::endl;
    // char* parameter_file "./lib/ForceField/Parameters.cif"
  cif_option_t compiler_options = cif_option_default();
  cexception_t inner;
  // std::cout << __FILE__ << std::endl;
  // CIF* parameters =
  //   new_cif_from_cif_file(parameter_file, compiler_options, &inner);
  // std::cout << cif_tag_index(parameters, "_rotag_atom_properties.type_symbol") << std::endl;
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
