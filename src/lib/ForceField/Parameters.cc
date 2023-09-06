#include "Parameters.h"

#include <cmath>
#include <iostream>
#include <memory>

#include <boost/filesystem.hpp>

extern "C" {
  #include "cif_compiler.h"
}

Parameters::Parameters(char* program_file_path) {
  boost::filesystem::path parameter_file =
      boost::filesystem::canonical(program_file_path).parent_path().parent_path() /
      boost::filesystem::path(__FILE__).parent_path() /
      "Parameters.cif";
  cif_option_t compiler_options = cif_option_default();
  cexception_t inner;

  CIF* parameters =
    new_cif_from_cif_file((char*) parameter_file.c_str(),
                          compiler_options,
                          &inner);
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
