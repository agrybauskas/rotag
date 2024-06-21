#ifndef SRC_LIB_CIFTAG_H_
#define SRC_LIB_CIFTAG_H_

#include <map>
#include <string>
#include <vector>

extern "C" {
  #include "datablock.h"
}

std::map<std::string, ssize_t> cif_tag_index_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags);

std::map<std::string, ssize_t> cif_value_length_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags);

#endif  // SRC_LIB_CIFTAG_H_
