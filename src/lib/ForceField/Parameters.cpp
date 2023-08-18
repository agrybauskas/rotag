#include "Parameters.h"

#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/join.hpp>

// #include "../../externals/codcif/"

Parameters::Parameters() {
  // NOTE: not sure if getenv() is ok here.
  std::vector<std::string> path = {
    std::getenv("ROTAG_SRC"), "src/lib/ForceField/Parameters.cif"
  };
  force_field_file = boost::join(path, "/");

  std::vector<std::string> data_identifier = {
    "_rotag_force_field.lj_k",
    "_rotag_force_field.c_k",
    "_rotag_force_field.h_k",
    "_rotag_force_field.t_k",
    "_rotag_force_field.cutoff_atom",
    "_rotag_force_field.cutoff_start",
    "_rotag_force_field.cutoff_end",
    "_rotag_atom_properties",
    "_rotag_lennard_jones",
    "_rotag_partial_charge",
    "_rotag_torsional_atom_names",
    "_rotag_torsional",
    "_rotag_h_bond",
    "_rotag_residue_atom_necessity",
    "_rotag_clear_hybridization",
    "_rotag_connectivity",
    "_rotag_hydrogen_names",
    "_rotag_symmetrical_atom_names",
    "_rotag_dihedral_angle",
    "_rotag_interaction_atom_names",
    "_rotag_mainchain_atom_names",
    "_rotag_sidechain_atom_names",
    "_rotag_rotatable_residue_names"
  };


}

void Parameters::_retrieve_constants() {

};

void Parameters::_retrieve_bond_data() {

};
