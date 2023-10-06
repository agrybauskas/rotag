#include "PDBxParser.h"

#include <iostream>

extern "C" {
  #include "cif_compiler.h"
}

#include "CIFTag.h"

AtomSite mmcif_to_atom_site(char* mmcif_file_path) {
  AtomSite atom_site = {};

  cif_option_t compiler_options = cif_option_default();
  cexception_t inner;
  CIF* mmcif = new_cif_from_cif_file(mmcif_file_path, compiler_options, &inner);

  const std::vector<std::string> cif_tags = {
    "_atom_site.group_pdb",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.cartn_x",
    "_atom_site.cartn_y",
    "_atom_site.cartn_z",
    "_atom_site.occupancy",
    "_atom_site.b_iso_or_equiv",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
    "_atom_site.pdbx_pdb_model_num",
    "_atom_site.rotag_selection_state",
    "_atom_site.rotag_selection_group"
  };

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(mmcif)) {
    std::map<std::string, ssize_t> cif_tag_index_lookup_table =
      cif_tag_index_lookup(datablock, cif_tags);
    std::map<std::string, ssize_t> cif_value_length_lookup_table =
      cif_value_length_lookup(datablock, cif_tags);

    // "_atom_site" category.
    for (int i = 0; i < cif_value_length_lookup_table["_atom_site.id"]; i++) {
      Atom atom = {};
      if(cif_tag_index_lookup_table["_atom_site.group_pdb"] > 0) {
        atom.group_pdb = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.group_pdb", i);
      }
      if(cif_tag_index_lookup_table["_atom_site.id"] > 0) {
        atom.id = cifvalue_to_long(
          datablock, cif_tag_index_lookup_table, "_atom_site.id", i);
      }
      if(cif_tag_index_lookup_table["_site.type_symbol"] > 0) {
        atom.type_symbol = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.type_symbol", i);
      }
      if(cif_tag_index_lookup_table[".label_atom_id"] > 0) {
        atom.label_atom_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.label_atom_id", i);
      }
      if(cif_tag_index_lookup_table[".label_alt_id"] > 0) {
        atom.label_alt_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.label_alt_id", i);
      }
      if(cif_tag_index_lookup_table[".label_comp_id"] > 0) {
        atom.label_comp_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.label_comp_id", i);
      }
      if(cif_tag_index_lookup_table[".label_asym_id"] > 0) {
        atom.label_asym_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.label_asym_id", i);
      }
      if(cif_tag_index_lookup_table[".label_entity_id"] > 0) {
        atom.label_entity_id = cifvalue_to_long(
          datablock,
          cif_tag_index_lookup_table,
          "_atom_site.label_entity_id",
          i);
      }
      if(cif_tag_index_lookup_table[".label_seq_id"] > 0) {
        atom.label_seq_id = cifvalue_to_long(
          datablock, cif_tag_index_lookup_table, "_atom_site.label_seq_id", i);
      }
      if(cif_tag_index_lookup_table["_site.cartn_x"] > 0) {
        atom.cartn_x = cifvalue_to_double(
          datablock, cif_tag_index_lookup_table, "_atom_site.cartn_x", i);
      }
      if(cif_tag_index_lookup_table["_site.cartn_y"] > 0) {
        atom.cartn_y = cifvalue_to_double(
          datablock, cif_tag_index_lookup_table, "_atom_site.cartn_y", i);
      }
      if(cif_tag_index_lookup_table["_site.cartn_z"] > 0) {
        atom.cartn_z = cifvalue_to_double(
          datablock, cif_tag_index_lookup_table, "_atom_site.cartn_z", i);
      }
      if(cif_tag_index_lookup_table["_atom_site.occupancy"] > 0) {
        atom.occupancy = cifvalue_to_double(
          datablock, cif_tag_index_lookup_table, "_atom_site.occupancy", i);
      }
      if(cif_tag_index_lookup_table["_iso_or_equiv"] > 0) {
        atom.b_iso_or_equiv = cifvalue_to_double(
          datablock,
          cif_tag_index_lookup_table,
          "_atom_site.b_iso_or_equiv",
          i);
      }
      if(cif_tag_index_lookup_table[".auth_seq_id"] > 0) {
        atom.auth_seq_id = cifvalue_to_long(
          datablock, cif_tag_index_lookup_table, "_atom_site.auth_seq_id", i);
      }
      if(cif_tag_index_lookup_table[".auth_comp_id"] > 0) {
        atom.auth_comp_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.auth_comp_id", i);
      }
      if(cif_tag_index_lookup_table[".auth_asym_id"] > 0) {
        atom.auth_asym_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.auth_asym_id", i);
      }
      if(cif_tag_index_lookup_table[".auth_atom_id"] > 0) {
        atom.auth_atom_id = cifvalue_to_string(
          datablock, cif_tag_index_lookup_table, "_atom_site.auth_atom_id", i);
      }
      if(cif_tag_index_lookup_table["_atom_site.pdbx_pdb_model_num"] > 0) {
        atom.pdbx_pdb_model_num = cifvalue_to_long(
          datablock,
          cif_tag_index_lookup_table,
          "_atom_site.pdbx_pdb_model_num",
          i);
      }
      if(cif_tag_index_lookup_table["_atom_site.rotag_selection_state"] > 0) {
        atom.selection_state = cifvalue_to_string(
          datablock,
          cif_tag_index_lookup_table,
          "_atom_site.rotag_selection_state",
          i);
      }
      if(cif_tag_index_lookup_table["_atom_site.rotag_selection_group"] > 0) {
        atom.selection_group = cifvalue_to_string(
          datablock,
          cif_tag_index_lookup_table,
          "_atom_site.rotag_selection_group",
          i);
      }

      atom_site[atom.id] = atom;
    }
  }

  return atom_site;
}

std::vector<unsigned long int>
  filter(AtomSite atom_site, Selector include, Selector exclude) {

  if (atom_site.empty()) {
    /* TODO: Error or warning. Message: no atom were loaded to the AtomSite data
       structure. */
  }

  return std::vector<unsigned long int>{};
}
