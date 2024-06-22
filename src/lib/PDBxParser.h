#ifndef SRC_LIB_PDBXPARSER_H_
#define SRC_LIB_PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>

extern "C" {
  #include "cif_compiler.h"
}

const std::vector<std::string> ATOM_SITE_TAGS = {
  // "_atom_site" category-related.
  "_atom_site.group_pdb",              // "ATOM" or "HETATM".
  "_atom_site.id",                     // Atom id.
  "_atom_site.type_symbol",            // Chemical element.
  "_atom_site.label_atom_id",          // Atom label.
  "_atom_site.label_alt_id",           // Related to alternative atom position.
  "_atom_site.label_comp_id",          // Residue name.
  "_atom_site.label_asym_id",          // Chain name.
  "_atom_site.label_entity_id",        // Molecular entity.
  "_atom_site.label_seq_id",           // Residue id.
  "_atom_site.cartn_x",                // Cartesian x coordinates of the atom.
  "_atom_site.cartn_y",                // Cartesian y coordinates of the atom.
  "_atom_site.cartn_z",                // Cartesian z coordinates of the atom.
  "_atom_site.occupancy",              // The fraction present in the site.
  "_atom_site.b_iso_or_equiv",         // Isotropic displacement.
  "_atom_site.auth_seq_id",            // Author's residue id.
  "_atom_site.auth_comp_id",           // Author's residue name.
  "_atom_site.auth_asym_id",           // Author's chain name.
  "_atom_site.auth_atom_id",           // Author's atom label.
  "_atom_site.pdbx_pdb_model_num",     // Model id.

  // Selection that is specific to rotag.
  "_atom_site.rotag_selection_state",  // Marks selection state: T, S ir H.
  "_atom_site.rotag_selection_group"   // Selection group id.
};

const std::vector<std::string> PARAMETERS_TAGS = {
  // "_rotag_parameters" category-related.
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

// typedef std::map<std::string, std::map<std::string, std::string>> AtomSite;
// typedef std::map<std::string, std::map<std::string, bool>> Selector;

// AtomSite mmcif_to_atom_site(char* mmcif_file_path);

// AtomSite
// filter(AtomSite atom_site, Selector include = {{}}, Selector exclude = {{}});

// std::vector<std::vector<std::string>>
// extract(AtomSite atom_site, std::vector<std::string> cif_tags);

// void mark_selection(
//   AtomSite* atom_site,
//   std::vector<std::string> target_atom_ids = {},
//   std::vector<std::string> selected_atom_ids = {});

#endif  // SRC_LIB_PDBXPARSER_H_
