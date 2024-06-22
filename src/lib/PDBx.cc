#include "PDBx.h"

PDBx::PDBx(CIF* cif, std::vector<std::string> select_tags) {
  std::vector<std::string> cif_tags = {};
  if(select_tags.size() > 0) {
    cif_tags = select_tags;
  } else {
    char** tags = datablock_tags(cif_datablock_list(cif));
    for(int i = 0; i < sizeof(tags); i++ ) {
      std::string cif_tag = tags[i];
      cif_tags.push_back(cif_tag);
    }
  }

  DATABLOCK* datablock;
  foreach_datablock(datablock, cif_datablock_list(cif)) {
  }
}

PDBx::~PDBx(){}
