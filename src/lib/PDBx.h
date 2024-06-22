#ifndef SRC_LIB_PDBX_H_
#define SRC_LIB_PDBX_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

extern "C" {
  #include "cif.h"
  #include "datablock.h"
}

#include "PDBxParser.h"

struct PDBXVALUE {
  std::string value_str;
  double value_float;
  int64_t value_int;

  operator std::string () const { return value_str; }
  operator double () const { return value_float; }
  operator int64_t () const { return value_int; }

  explicit PDBXVALUE(CIFVALUE* x) {

    this->value_str = "A";
    this->value_float = 1.0;
    this->value_int = 2;
  }

  ~PDBXVALUE() {};
};

class PDBx {
 private:
  std::map<std::string, std::vector<PDBXVALUE>> data;
  std::map<std::string, std::vector<std::string>> categories;
  std::map<std::string, std::vector<double>> category_order;
  std::map<std::string, bool> loops;

 public:
  explicit PDBx(CIF* cif, std::vector<std::string> select_tags = {});
  ~PDBx();
  void values(std::string cif_tag);
};

#endif  // SRC_LIB_PDBX_H_
