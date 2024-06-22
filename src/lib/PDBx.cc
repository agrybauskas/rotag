#include "PDBx.h"

PDBx::PDBx(CIF* cif, std::vector<std::string> select_tags) {
  std::vector<std::string> cif_tags = {};
  if (select_tags.size() > 0) {
    cif_tags = select_tags;
  } else {
    char** tags = datablock_tags(cif_datablock_list(cif));
    for (int i = 0; i < (int) sizeof(tags); i++) {
      std::string cif_tag = tags[i];
      cif_tags.push_back(cif_tag);
    }
  }

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(cif)) {
    const ssize_t* cif_value_lengths = datablock_value_lengths(datablock);
    for (const std::string &cif_tag : cif_tags) {
      const ssize_t cif_tag_index =
        datablock_tag_index(datablock, (char*) cif_tag.c_str());

      for (ssize_t i = 0; i < cif_value_lengths[cif_tag_index]; i++) {
        datablock_cifvalue(datablock, cif_tag_index, i);
      }
    }
  }
}

PDBx::~PDBx() {}

void PDBx::values(std::string cif_tag) {}
