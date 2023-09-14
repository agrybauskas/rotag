#include "Parameters.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <boost/filesystem.hpp>

extern "C" {
  #include "cif_compiler.h"
  #include "datablock.h"
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

  const std::vector<std::string> atom_properties_items = {
    "_rotag_atom_properties.type_symbol",
    "_rotag_atom_properties.hybridization",
    "_rotag_atom_properties.covalent_radius_value",
    "_rotag_atom_properties.covalent_radius_error",
    "_rotag_atom_properties.vdw_radius",
    "_rotag_atom_properties.lone_pair_count",
    "_rotag_atom_properties.valence"
  };

  for(const std::string &atom_properties_item: atom_properties_items) {
    DATABLOCK* datablock;
    foreach_datablock(datablock, cif_datablock_list(parameters)) {
      const ssize_t tag_index =
        datablock_tag_index(datablock, (char*) atom_properties_item.c_str());
    }
  }
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
