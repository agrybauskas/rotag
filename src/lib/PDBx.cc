#include "PDBx.h"

PDBx::PDBx(CIF* cif, std::vector<std::string> tags) {
  this->data = cif;
}

PDBx::~PDBx(){}
