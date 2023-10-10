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

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(mmcif)) {
    std::map<std::string, ssize_t> cif_tag_index_lookup_table =
      cif_tag_index_lookup(datablock, ATOM_SITE_TAGS);
    std::map<std::string, ssize_t> cif_value_length_lookup_table =
      cif_value_length_lookup(datablock, ATOM_SITE_TAGS);

    for (const std::string &cif_tag : ATOM_SITE_TAGS) {
      for (int i = 0; i < cif_value_length_lookup_table["_atom_site.id"]; i++) {
        if(cif_tag_index_lookup_table[cif_tag] > 0) {
          atom_site["_atom_site.id"][cif_tag] = value_scalar(datablock_cifvalue(
            datablock, cif_tag_index_lookup_table[cif_tag], i));
        }
      }
    }
  }

  return atom_site;
}

std::vector<std::string>
  filter(AtomSite atom_site, Selector include, Selector exclude) {

  if (atom_site.empty()) {
    // TODO: Error or warning. Message: no atoms were loaded from "_atom_site".
  }

  std::vector<std::string> filtered_atom_ids = {};
  for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
    bool keep_atom = true;
    // for () {

    // }
  }

  return filtered_atom_ids;
}
