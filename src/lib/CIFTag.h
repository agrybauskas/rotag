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

double cifvalue_to_double(
  DATABLOCK* datablock,
  std::map<std::string, ssize_t> cif_tag_index_lookup_table,
  std::string cif_tag,
  size_t index = 0);

std::string cifvalue_to_string(
  DATABLOCK* datablock,
  std::map<std::string, ssize_t> cif_tag_index_lookup_table,
  std::string cif_tag,
  size_t index = 0);

#endif  // SRC_LIB_CIFTAG_H_
