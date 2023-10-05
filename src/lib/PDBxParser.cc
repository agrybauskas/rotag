#include "PDBxParser.h"

AtomSite mmcif_to_atom_site(char* mmcif_file_path) {
  AtomSite atom_site = {};

  cif_option_t compiler_options = cif_option_default();
  cexception_t inner;
  CIF* mmcif = new_cif_from_cif_file(mmcif_file_path, compiler_options, &inner);

  return atom_site;
}
