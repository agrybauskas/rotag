#ifndef SRC_LIB_CIFTAG_H_
#define SRC_LIB_CIFTAG_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

extern "C" {
  #include "cifvalue.h"
  #include "ciflist.h"
  #include "ciftable.h"
  #include "datablock.h"
}

// Re-defining CIFVALUE in order to have proper explicit conversions.
struct CIFVALUE {
  union {
    char *str;
    struct CIFLIST *l;
    struct CIFTABLE *t;
  } v;
  cif_value_type_t type;
};

std::map<std::string, ssize_t> cif_tag_index_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags);

std::map<std::string, ssize_t> cif_value_length_lookup(
  DATABLOCK* datablock, std::vector<std::string> cif_tags);

#endif  // SRC_LIB_CIFTAG_H_
