#include "CIFTag.h"

#include <iostream>

std::map<std::string, ssize_t> cif_tag_index_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags) {
  std::map<std::string, ssize_t> cif_tag_index_lookup;
  for (const std::string &cif_tag : cif_tags) {
    cif_tag_index_lookup[cif_tag] =
      datablock_tag_index(datablock, (char*) cif_tag.c_str());
  }
  return cif_tag_index_lookup;
}

std::map<std::string, ssize_t> cif_value_length_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags) {
  std::map<std::string, ssize_t> cif_value_length_lookup;
  for (const std::string &cif_tag : cif_tags) {
    cif_value_length_lookup[cif_tag] =
      datablock_value_lengths(datablock)
        [datablock_tag_index(datablock, (char*) cif_tag.c_str())];
  }
  return cif_value_length_lookup;
}
