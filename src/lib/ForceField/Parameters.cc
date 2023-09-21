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

  const std::vector<std::string> cif_tags = {
    "_rotag_force_field.lj_k",
    "_rotag_force_field.c_k",
    "_rotag_force_field.h_k",
    "_rotag_force_field.t_k",
    "_rotag_force_field.cutoff_atom",
    "_rotag_force_field.cutoff_start",
    "_rotag_force_field.cutoff_end",
    "_rotag_atom_properties.type_symbol",
    "_rotag_atom_properties.hybridization",
    "_rotag_atom_properties.covalent_radius_value",
    "_rotag_atom_properties.covalent_radius_error",
    "_rotag_atom_properties.vdw_radius",
    "_rotag_atom_properties.lone_pair_count",
    "_rotag_atom_properties.valence"
  };

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(parameters)) {
    // std::map<std::string, std::map<unsigned long long int, char>> cif_data_lookup_table;
    for(const std::string &cif_tag: cif_tags) {
      const ssize_t tag_index =
        datablock_tag_index(datablock, (char*) cif_tag.c_str());
      const ssize_t tag_value_lengths =
        datablock_value_lengths(datablock)[tag_index];
      for(int i = 0; i < tag_value_lengths; i++) {

        // "_rotag_force_field" category.
        if(cif_tag == "_rotag_force_field.lj_k") {
          this->lj_k =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.c_k") {
          this->c_k =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.h_k") {
          this->h_k =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.t_k") {
          this->t_k =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.cutoff_atom") {
          this->cutoff_atom =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.cutoff_start") {
          this->cutoff_start =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));
        } else if(cif_tag == "_rotag_force_field.cutoff_end") {
          this->cutoff_end =
            atof(value_scalar(datablock_cifvalue(datablock, tag_index, i)));

        // "_rotag_atom_properties" category.
        } else if(cif_tag == "_rotag_atom_properties.type_symbol") {
        }
      }
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
