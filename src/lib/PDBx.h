#ifndef SRC_LIB_PDBX_H_
#define SRC_LIB_PDBX_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>

extern "C" {
    #include "cif_compiler.h"
    #include "cifvalue.h"
    #include "datablock.h"
}

struct PDBXVALUE {
    std::string value_str;
    double value_float;
    int64_t value_int;

    operator std::string () const { return value_str; }
    operator double () const { return value_float; }
    operator int64_t () const { return value_int; }

    explicit PDBXVALUE(CIFVALUE* cif_value) {
        cif_value_type_t type = value_type(cif_value);
        switch (type) {
            case CIF_INT:
                this->value_int = atoi(value_scalar(cif_value));
                this->value_str = value_scalar(cif_value);
                break;
            case CIF_FLOAT:
                this->value_float = atof(value_scalar(cif_value));
                this->value_str = value_scalar(cif_value);
                break;
            case CIF_UQSTRING:
            case CIF_SQSTRING:
            case CIF_DQSTRING:
            case CIF_SQ3STRING:
            case CIF_DQ3STRING:
            case CIF_TEXT:
                this->value_str = value_scalar(cif_value);
                break;
            // HACK: look at these cases: CIF_UNKNOWN, CIF_NON_EXISTANT,
            // CIF_LIST, CIF_TABLE, last_CIF_VALUE.
            default:
                this->value_str = "";
                break;
        }
    }

    ~PDBXVALUE() {};
};

class PDBx {
 private:
    std::map<std::string, std::vector<PDBXVALUE>> data;
    std::map<std::string, int64_t> cif_tag_order;
    std::map<std::string, int64_t> in_loop;

 public:
    explicit PDBx(CIF* cif, std::vector<std::string> select_tags = {});
    ~PDBx();
    std::vector<PDBXVALUE> values(std::string cif_tag);
};

#endif  // SRC_LIB_PDBX_H_
