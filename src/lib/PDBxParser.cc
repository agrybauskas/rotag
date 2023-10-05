#include "PDBxParser.h"

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
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
    "_atom_site.pdbx_PDB_model_num",
    "_atom_site.rotag_selection_state",
    "_atom_site.rotag_selection_group"
  };

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(mmcif)) {
    std::map<std::string, ssize_t> cif_tag_index_lookup_table =
      cif_tag_index_lookup(datablock, cif_tags);
    std::map<std::string, ssize_t> cif_value_length_lookup_table =
      cif_value_length_lookup(datablock, cif_tags);
  }

  return atom_site;
}
