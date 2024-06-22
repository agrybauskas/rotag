#ifndef SRC_LIB_PDBX_H_
#define SRC_LIB_PDBX_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

extern "C" {
  #include "cif.h"
  #include "cif_compiler.h"
  #include "datablock.h"
}

#include "PDBxParser.h"

class PDBx {
  private:
    CIF* data;
    std::vector<double> order;
    std::map<std::string, std::vector<std::string>> group;
    std::map<std::string, std::string> id;

  public:
    PDBx(CIF* cif, std::vector<std::string> select_tags);
    ~PDBx();
};

#endif  // SRC_LIB_PDBX_H_
