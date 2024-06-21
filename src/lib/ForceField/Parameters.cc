#include "Parameters.h"

#include <cmath>
#include <memory>

#include <boost/filesystem.hpp>

extern "C" {
  #include "cif_compiler.h"
}

#include "../CIFTag.h"
#include "../Combinatorics.h"

Parameters::Parameters(char* program_file_path) {
  boost::filesystem::path parameter_file =
    boost::filesystem::canonical(program_file_path).parent_path().parent_path()/
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
    "_rotag_atom_properties.valence",
    "_rotag_lennard_jones.type_symbol_1",
    "_rotag_lennard_jones.type_symbol_2",
    "_rotag_lennard_jones.sigma",
    "_rotag_lennard_jones.epsilon",
    "_rotag_partial_charge.label_comp_id",
    "_rotag_partial_charge.label_atom_id",
    "_rotag_partial_charge.value",
    "_rotag_torsional_atom_names.label_comp_id",
    "_rotag_torsional_atom_names.label_atom_id",
    "_rotag_torsional_atom_names.alt_atom_name",
    "_rotag_torsional.label_atom_1_id",
    "_rotag_torsional.label_atom_2_id",
    "_rotag_torsional.label_atom_3_id",
    "_rotag_torsional.label_atom_4_id",
    "_rotag_torsional.epsilon",
    "_rotag_torsional.phase",
    "_rotag_torsional.gamma",
    "_rotag_h_bond.type_symbol",
    "_rotag_h_bond.sigma",
    "_rotag_h_bond.epsilon",
    "_rotag_residue_atom_necessity.label_comp_id",
    "_rotag_residue_atom_necessity.label_atom_id",
    "_rotag_residue_atom_necessity.value",
    "_rotag_clear_hybridization.label_comp_id",
    "_rotag_clear_hybridization.label_atom_id",
    "_rotag_clear_hybridization.type",
    "_rotag_connectivity.label_comp_id",
    "_rotag_connectivity.label_atom_1_id",
    "_rotag_connectivity.label_atom_2_id",
    "_rotag_hydrogen_names.label_comp_id",
    "_rotag_hydrogen_names.label_atom_id",
    "_rotag_hydrogen_names.label_hydrogen_atom_id",
    "_rotag_symmetrical_atom_names.label_comp_id",
    "_rotag_symmetrical_atom_names.label_atom_1_id",
    "_rotag_symmetrical_atom_names.label_atom_2_id",
    "_rotag_dihedral_angle.label_comp_id",
    "_rotag_dihedral_angle.angle",
    "_rotag_dihedral_angle.range_from",
    "_rotag_dihedral_angle.range_to",
    "_rotag_dihedral_angle.step",
    "_rotag_dihedral_angle.type",
    "_rotag_interaction_atom_names.label_atom_id",
    "_rotag_mainchain_atom_names.label_atom_id",
    "_rotag_sidechain_atom_names.label_atom_id",
    "_rotag_rotatable_residue_names.label_comp_id"
  };

  // DATABLOCK* datablock;
  // foreach_datablock(datablock, cif_datablock_list(parameters)) {
  //   std::map<std::string, ssize_t> cif_tag_index_lookup_table;
  //   std::map<std::string, ssize_t> cif_value_length_lookup_table;

  //   // Generating lookup tables first.
  //   for (const std::string &cif_tag : cif_tags) {
  //     cif_tag_index_lookup_table[cif_tag] =
  //       datablock_tag_index(datablock, (char*) cif_tag.c_str());
  //     cif_value_length_lookup_table[cif_tag] =
  //       datablock_value_lengths(datablock)[cif_tag_index_lookup_table[cif_tag]];
  //   }

  //   // Parsing tags per case basis.
  //   // "_rotag_force_field" category.
  //   this->lj_k = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.lj_k");
  //   this->c_k = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.c_k");
  //   this->h_k = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.h_k");
  //   this->t_k = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.t_k");
  //   this->cutoff_atom = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.cutoff_atom");
  //   this->cutoff_start = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.cutoff_start");
  //   this->cutoff_end = cifvalue_to_double(
  //     datablock, cif_tag_index_lookup_table, "_rotag_force_field.cutoff_end");

  //   /* NOTE: "codcif" parser should catch errors if the length of tag values
  //      does not have the same size in the loop so, it is enough to choose any
  //      column for iterating. */

  //   // "_rotag_atom_properties" category.
  //   for (int i = 0;
  //        i <cif_value_length_lookup_table["_rotag_atom_properties.type_symbol"];
  //        i++) {
  //     std::string type_symbol = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_atom_properties.type_symbol",
  //       i);
  //     std::string hybridization = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_atom_properties.hybridization",
  //       i);

  //     this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value =
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_atom_properties.covalent_radius_value",
  //         i);
  //     this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].error =
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_atom_properties.covalent_radius_error",
  //         i);
  //     this->ATOM_PROPERTIES[type_symbol].vdw_radius =
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_atom_properties.vdw_radius",
  //         i);
  //     this->ATOM_PROPERTIES[type_symbol].lone_pair_count =
  //       std::stoi(cifvalue_to_string(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_atom_properties.lone_pair_count",
  //         i));
  //     this->ATOM_PROPERTIES[type_symbol].valence =
  //       std::stoi(cifvalue_to_string(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_atom_properties.valence",
  //         i));

  //     // Used for edge length of the grid during cubing procedure algorithm.
  //     if (this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value >
  //         this->max_connection_length) {
  //       this->max_connection_length =
  //         this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value;
  //     }

  //     /* Covalent radii are stored in list as it will be much more easier and
  //        cleaner to calculate bond length combinations. */
  //     // NOTE: values and errors should have the count in list.
  //     if (this->COVALENT_RADII_VALUES.count(type_symbol) > 0) {
  //       this->COVALENT_RADII_VALUES[type_symbol].push_back(
  //         this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value);
  //       this->COVALENT_RADII_ERRORS[type_symbol].push_back(
  //         this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].error);
  //     } else {
  //       this->COVALENT_RADII_VALUES[type_symbol] = std::vector<double>{
  //         this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].value
  //       };
  //       this->COVALENT_RADII_ERRORS[type_symbol] = std::vector<double>{
  //         this->ATOM_PROPERTIES[type_symbol].covalent_radius[hybridization].error
  //       };
  //     }
  //   }

  //   // "_rotag_partial_charge" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_lennard_jones.type_symbol_1"];
  //        i++) {
  //     std::string type_symbol_1 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_lennard_jones.type_symbol_1",
  //       i);
  //     std::string type_symbol_2 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_lennard_jones.type_symbol_2",
  //       i);
  //     double sigma = cifvalue_to_double(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_lennard_jones.sigma",
  //       i);
  //     double epsilon = cifvalue_to_double(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_lennard_jones.epsilon",
  //       i);

  //     this->LENNARD_JONES[type_symbol_1][type_symbol_2] = {sigma, epsilon};
  //     this->LENNARD_JONES[type_symbol_2][type_symbol_1] = {sigma, epsilon};
  //   }

  //   // "_rotag_partial_charge" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_partial_charge.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_partial_charge.label_comp_id",
  //       i);
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_partial_charge.label_atom_id",
  //       i);
  //     double partial_charge_value = cifvalue_to_double(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_partial_charge.value",
  //       i);

  //     this->PARTIAL_CHARGE[residue_name][atom_name].value =
  //       partial_charge_value;
  //   }

  //   // "_rotag_torsional_atom_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_torsional_atom_names.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional_atom_names.label_comp_id",
  //       i);
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional_atom_names.label_atom_id",
  //       i);
  //     std::string alt_atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional_atom_names.alt_atom_name",
  //       i);

  //     this->TORSIONAL_ATOM_NAMES[residue_name][atom_name].alt_name =
  //       alt_atom_name;
  //   }

  //   // "_rotag_torsional" category,
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_torsional.label_atom_1_id"];
  //        i++) {
  //     std::string atom_name_1 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional.label_atom_1_id",
  //       i);
  //     std::string atom_name_2 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional.label_atom_2_id",
  //       i);
  //     std::string atom_name_3 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional.label_atom_3_id",
  //       i);
  //     std::string atom_name_4 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_torsional.label_atom_4_id",
  //       i);

  //     std::string key_forward =
  //       atom_name_1 + "," + atom_name_2 + "," + atom_name_3 + "," + atom_name_4;
  //     std::string key_reverse =
  //       atom_name_4 + "," + atom_name_3 + "," + atom_name_2 + "," + atom_name_1;

  //     double epsilon = cifvalue_to_double(
  //       datablock, cif_tag_index_lookup_table, "_rotag_torsional.epsilon", i);
  //     double phase = cifvalue_to_double(
  //       datablock, cif_tag_index_lookup_table, "_rotag_torsional.phase", i);
  //     double gamma = cifvalue_to_double(
  //       datablock, cif_tag_index_lookup_table, "_rotag_torsional.gamma", i);

  //     this->TORSIONAL[key_forward] = {epsilon, phase, gamma};
  //     this->TORSIONAL[key_reverse] = {epsilon, phase, gamma};
  //   }

  //   // "_rotag_h_bond" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_h_bond.type_symbol"];
  //        i++) {
  //     std::string type_symbol = cifvalue_to_string(
  //       datablock, cif_tag_index_lookup_table, "_rotag_h_bond.type_symbol", i);
  //     double sigma = cifvalue_to_double(
  //       datablock, cif_tag_index_lookup_table, "_rotag_h_bond.sigma", i);
  //     double epsilon = cifvalue_to_double(
  //       datablock, cif_tag_index_lookup_table, "_rotag_h_bond.epsilon", i);

  //     this->H_BOND[type_symbol] = {sigma, epsilon};
  //   }

  //   // "_rotag_residue_atom_necessity" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_residue_atom_necessity.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_residue_atom_necessity.label_comp_id",
  //       i);
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_residue_atom_necessity.label_atom_id",
  //       i);
  //     std::string necessity_value = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_residue_atom_necessity.value",
  //       i);

  //     this->RESIDUE_ATOM_NECESSITY[residue_name][atom_name] =
  //       necessity_value == "mandatory" ? true : false;
  //   }

  //   // "_rotag_clear_hybridization" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_clear_hybridization.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_clear_hybridization.label_comp_id",
  //       i);
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_clear_hybridization.label_atom_id",
  //       i);
  //     std::string type = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_clear_hybridization.type",
  //       i);

  //     this->CLEAR_HYBRIDIZATION[residue_name][atom_name].type = type;
  //   }

  //   // "_rotag_connectivity" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_connectivity.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_connectivity.label_comp_id",
  //       i);
  //     std::string atom_name_1 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_connectivity.label_atom_1_id",
  //       i);
  //     std::string atom_name_2 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_connectivity.label_atom_2_id",
  //       i);

  //     this->CONNECTIVITY[residue_name][atom_name_1].push_back(atom_name_2);
  //   }

  //   // "_rotag_hydrogen_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_hydrogen_names.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_hydrogen_names.label_comp_id",
  //       i);
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_hydrogen_names.label_atom_id",
  //       i);
  //     std::string hydrogen_atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_hydrogen_names.label_hydrogen_atom_id",
  //       i);

  //     this->HYDROGEN_NAMES[residue_name][atom_name].push_back(
  //       hydrogen_atom_name);
  //   }

  //   // "_rotag_symmetrical_atom_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_symmetrical_atom_names.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_symmetrical_atom_names.label_comp_id",
  //       i);
  //     std::string atom_name_1 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_symmetrical_atom_names.label_atom_1_id",
  //       i);
  //     std::string atom_name_2 = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_symmetrical_atom_names.label_atom_2_id",
  //       i);

  //     this->SYMMETRICAL_ATOM_NAMES[residue_name][atom_name_1].push_back(
  //       atom_name_2);
  //     this->SYMMETRICAL_ATOM_NAMES[residue_name][atom_name_2].push_back(
  //       atom_name_1);
  //   }

  //   // "_rotag_dihedral_angle" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_dihedral_angle.label_comp_id"];
  //        i++) {
  //     std::string residue_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_dihedral_angle.label_comp_id",
  //       i);
  //     std::string angle = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_dihedral_angle.angle",
  //       i);
  //     double range_from =
  //       (std::string) cifvalue_to_string(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.range_from",
  //         i) == "." ?
  //       NULL :
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.range_from",
  //         i);
  //     double range_to =
  //       (std::string) cifvalue_to_string(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.range_to",
  //         i) == "." ?
  //       NULL :
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.range_to",
  //         i);
  //     double step =
  //       (std::string) cifvalue_to_string(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.step",
  //         i) == "." ?
  //       NULL :
  //       cifvalue_to_double(
  //         datablock,
  //         cif_tag_index_lookup_table,
  //         "_rotag_dihedral_angle.step",
  //         i);
  //     std::string type = cifvalue_to_string(
  //       datablock, cif_tag_index_lookup_table, "_rotag_dihedral_angle.type", i);

  //     this->DIHEDRAL_ANGLE[residue_name][angle] = {
  //       range_from,
  //       range_to,
  //       step,
  //       type
  //     };
  //   }

  //   // "_rotag_interaction_atom_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_interaction_atom_names.label_atom_id"];
  //        i++) {
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_interaction_atom_names.label_atom_id",
  //       i);

  //     this->INTERACTION_ATOM_NAMES.push_back(atom_name);
  //   }

  //   // "_rotag_mainchain_atom_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_mainchain_atom_names.label_atom_id"];
  //        i++) {
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_mainchain_atom_names.label_atom_id",
  //       i);

  //     this->MAINCHAIN_ATOM_NAMES.push_back(atom_name);
  //   }

  //   // "_rotag_sidechain_atom_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_sidechain_atom_names.label_atom_id"];
  //        i++) {
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_sidechain_atom_names.label_atom_id",
  //       i);

  //     this->SIDECHAIN_ATOM_NAMES.push_back(atom_name);
  //   }

  //   // "_rotag_rotatable_residue_names" category.
  //   for (int i = 0;
  //        i < cif_value_length_lookup_table["_rotag_rotatable_residue_names.label_comp_id"];
  //        i++) {
  //     std::string atom_name = cifvalue_to_string(
  //       datablock,
  //       cif_tag_index_lookup_table,
  //       "_rotag_rotatable_residue_names.label_comp_id",
  //       i);

  //     this->ROTATABLE_RESIDUE_NAMES.push_back(atom_name);
  //   }
  // }

  // // Arginine model is used for calculating the interaction cutoff.
  // this->max_interaction_length =
  //   9 * this->ATOM_PROPERTIES["C"].covalent_radius["sp3"].value +
  //   2 * this->ATOM_PROPERTIES["N"].covalent_radius["sp3"].value +
  //       this->ATOM_PROPERTIES["H"].covalent_radius["sp3"].value;

  // // Determining min/max bond lengths for specific bond types.
  // // Single bonds.
  // const std::vector<std::string> atom_symbols_single ={"C", "H", "N", "O", "S"};
  // for (const std::string &first_atom_symbol : atom_symbols_single) {
  //   for (const std::string &second_atom_symbol : atom_symbols_single) {
  //     this->BOND_TYPE["single"][first_atom_symbol][second_atom_symbol] = {
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp3"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp3"].value -
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp3"].error -
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp3"].error,
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp3"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp3"].value +
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp3"].error +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp3"].error
  //     };
  //   }
  // }

  // // Double bonds.
  // const std::vector<std::string> atom_symbols_double = {"C", "N", "O"};
  // for (const std::string &first_atom_symbol : atom_symbols_double) {
  //   for (const std::string &second_atom_symbol : atom_symbols_double) {
  //     this->BOND_TYPE["double"][first_atom_symbol][second_atom_symbol] = {
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp2"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp2"].value -
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp2"].error -
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp2"].error,
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp2"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp2"].value +
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp2"].error +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp2"].error
  //     };
  //   }
  // }

  // // Triple bonds.
  // const std::vector<std::string> atom_symbols_triple = {"C"};
  // for (const std::string &first_atom_symbol : atom_symbols_triple) {
  //   for (const std::string &second_atom_symbol : atom_symbols_triple) {
  //     this->BOND_TYPE["triple"][first_atom_symbol][second_atom_symbol] = {
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp"].value -
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp"].error -
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp"].error,
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp"].value +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp"].value +
  //       this->ATOM_PROPERTIES[first_atom_symbol].covalent_radius["sp"].error +
  //       this->ATOM_PROPERTIES[second_atom_symbol].covalent_radius["sp"].error
  //     };
  //   }
  // }

  // // Precalculating covalent bond combinations.
  // std::vector<std::string> atom_symbols = {};
  // for (std::map<std::string, AtomProperties>::iterator it_i =
  //        this->ATOM_PROPERTIES.begin();
  //      it_i != this->ATOM_PROPERTIES.end();
  //      ++it_i) {
  //   std::string first_atom_symbol = it_i->first;
  //   for (std::map<std::string, AtomProperties>::iterator it_j =
  //          this->ATOM_PROPERTIES.begin();
  //        it_j != this->ATOM_PROPERTIES.end();
  //        ++it_j) {
  //     std::string second_atom_symbol = it_j->first;
  //     std::vector<std::vector<double>> length_combinations = {};
  //     std::vector<std::vector<double>> error_combinations = {};

  //     permutation(
  //       2,
  //       std::vector<std::vector<double>>({
  //         this->COVALENT_RADII_VALUES[first_atom_symbol],
  //         this->COVALENT_RADII_VALUES[second_atom_symbol],
  //       }),
  //       &length_combinations);
  //     permutation(2,
  //         std::vector<std::vector<double>>({
  //         this->COVALENT_RADII_ERRORS[first_atom_symbol],
  //         this->COVALENT_RADII_ERRORS[second_atom_symbol],
  //             }),
  //       &error_combinations);

  //     this->COVALENT_BOND_COMBINATIONS[first_atom_symbol][second_atom_symbol].values =
  //       length_combinations;
  //     this->COVALENT_BOND_COMBINATIONS[first_atom_symbol][second_atom_symbol].errors =
  //       error_combinations;
  //   }
  // }
}

Parameters::~Parameters() {}

double Parameters::epsilon() {
  double epsilon = 1.0;
  while ((1.0 + 0.5 * epsilon) != 1.0) {
    epsilon = 0.5 * epsilon;
  }
  return epsilon;
}

double Parameters::pi() {
  return 4 * std::atan2(1, 1);
}
