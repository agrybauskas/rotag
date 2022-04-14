#include "Parameters.h"

#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/join.hpp>

#include "../PDBxParser.h"

Parameters::Parameters() {
  // TODO: have to implement relative path to current module.
  force_field_file = "src/lib/ForceField/Parameters.cif";

  std::vector<std::string> data_identifier = {
    "_rotag_force_field.lj_k","_rotag_force_field.c_k", "_rotag_force_field.h_k",
    "_rotag_force_field.t_k", "_rotag_force_field.cutoff_atom",
    "_rotag_force_field.cutoff_start", "_rotag_force_field.cutoff_end",
    "_rotag_atom_properties", "_rotag_lennard_jones", "_rotag_partial_charge",
    "_rotag_partial_charge", "_rotag_torsional_atom_names", "_rotag_torsional",
    "_rotag_h_bond","_rotag_residue_atom_necessity","_rotag_clear_hybridization",
    "_rotag_connectivity", "_rotag_hydrogen_names",
    "_rotag_symmetrical_atom_names", "_rotag_dihedral_angle",
    "_rotag_interaction_atom_names", "_rotag_mainchain_atom_names",
    "_rotag_sidechain_atom_names", "_rotag_rotatable_residue_names"
  };

  obtain_pdbx_data(force_field_file, data_identifier);
}

void Parameters::_retrieve_constants() {

};

void Parameters::_retrieve_bond_data() {

};
