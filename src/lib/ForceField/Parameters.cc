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
    std::map<std::string, ssize_t> cif_tag_index_lookup_table;
    std::map<std::string, ssize_t> cif_value_length_lookup_table;

    // Generating lookup tables first.
    for(const std::string &cif_tag: cif_tags) {
      cif_tag_index_lookup_table[cif_tag] =
        datablock_tag_index(datablock, (char*) cif_tag.c_str());
      cif_value_length_lookup_table[cif_tag] =
        datablock_value_lengths(datablock)[cif_tag_index_lookup_table[cif_tag]];
    }

    // Parsing tags per case basis.
    // "_rotag_force_field" category.
    this->lj_k =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.lj_k"], 0)));
    this->c_k =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.c_k"], 0)));
    this->h_k =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.h_k"], 0)));
    this->t_k =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.t_k"], 0)));
    this->cutoff_atom =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.cutoff_atom"], 0)));
    this->cutoff_start =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.cutoff_start"], 0)));
    this->cutoff_end =
      atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_force_field.cutoff_end"], 0)));

    /* NOTE: "codcif" parser should catch errors if the length of tag values does
       not have the same size in the loop so, it is enough to choose any column
       for iterating. */

    // "_rotag_atom_properties" category.
    for(int i = 0; i < cif_value_length_lookup_table["_rotag_atom_properties.type_symbol"]; i++ ) {
      std::string type_symbol =
        value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.type_symbol"], i));
      std::string hybridization =
        value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.hybridization"], i));

      this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.covalent_radius_value"], i)));
      this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].error =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.covalent_radius_error"], i)));
      this->ATOM_PROPERTIES[type_symbol].vdw_radius =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.vdw_radius"], i)));
      this->ATOM_PROPERTIES[type_symbol].lone_pair_count =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.lone_pair_count"], i)));
      this->ATOM_PROPERTIES[type_symbol].valence =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_atom_properties.valence"], i)));

      // Used for edge length of the grid during cubing procedure algorithm.
      if(this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value >
         this->max_connection_length) {
        this->max_connection_length =
          this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value;
      }
    }

    // "_rotag_partial_charge" category.
    for(int i = 0; i < cif_value_length_lookup_table["_rotag_lennard_jones.type_symbol_1"]; i++ ) {
      std::string type_symbol_1 =
        value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_lennard_jones.type_symbol_1"], i));
      std::string type_symbol_2 =
        value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_lennard_jones.type_symbol_2"], i));
      double sigma =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_lennard_jones.sigma"], i)));
      double epsilon =
        atof(value_scalar(datablock_cifvalue(datablock, cif_tag_index_lookup_table["_rotag_lennard_jones.epsilon"], i)));

      this->LENNARD_JONES[type_symbol_1][type_symbol_2].sigma = sigma;
      this->LENNARD_JONES[type_symbol_2][type_symbol_1].sigma = sigma;
      this->LENNARD_JONES[type_symbol_1][type_symbol_2].epsilon = epsilon;
      this->LENNARD_JONES[type_symbol_2][type_symbol_1].epsilon = epsilon;
    }

    // "_rotag_partial_charge" category.

    // "_rotag_partial_charge" category.

    /* "_rotag_torsional_atom_names" and "_rotag_torsional" categories. Combining
       them to more useful structure. */

    // "_rotag_h_bond" category.

    // "_rotag_residue_atom_necessity" category.

    // "_rotag_clear_hybridization" category.

    // "_rotag_connectivity" category.

    // "_rotag_hydrogen_names" category.

    // "_rotag_symmetrical_atom_names" category.

    // "_rotag_dihedral_angle" category.

    // "_rotag_interaction_atom_names" category.

    // "_rotag_mainchain_atom_names" category.

    // "_rotag_sidechain_atom_names" category.

    // "_rotag_rotatable_residue_names" category.
  }

  // Arginine model is used for calculating the interaction cutoff.
  this->max_interaction_length =
      9 * this->ATOM_PROPERTIES["C"].covalent_radius["sp3"].value +
      2 * this->ATOM_PROPERTIES["N"].covalent_radius["sp3"].value +
          this->ATOM_PROPERTIES["H"].covalent_radius["sp3"].value;
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
