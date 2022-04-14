#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>

class Parameters {
  private:
    std::string force_field_file;

    void _retrieve_constants();
    void _retrieve_atom_data();
    void _retrieve_bond_data();
    void _retrieve_force_field_data();

  public:
    Parameters();
};

#endif
