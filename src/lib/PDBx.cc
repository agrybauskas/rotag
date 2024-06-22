#include "PDBx.h"

PDBx::PDBx(CIF* cif, std::vector<std::string> select_tags) {
    std::vector<std::string> cif_tags = {};
    if (select_tags.size() > 0) {
        cif_tags = select_tags;
    } else {
        char** tags = datablock_tags(cif_datablock_list(cif));
        for (int i = 0; i < static_cast<int>(sizeof(tags)); i++) {
            std::string cif_tag = tags[i];
            cif_tags.push_back(cif_tag);
        }
    }

    DATABLOCK* datablock;
    foreach_datablock(datablock, cif_datablock_list(cif)) {
        const ssize_t* cif_value_lengths = datablock_value_lengths(datablock);
        const int* are_in_loops = datablock_in_loop(datablock);
        for (const std::string &cif_tag : cif_tags) {
            const ssize_t cif_tag_index =
                datablock_tag_index(datablock,
                                    const_cast<char*>(cif_tag.c_str()));
            const int is_in_loop = are_in_loops[cif_tag_index];
            std::cout << cif_tag << std::endl;
            std::cout << "    is_in_loop: " << is_in_loop << std::endl;
            for (ssize_t i = 0; i < cif_value_lengths[cif_tag_index]; i++) {
                // PDBXVALUE pdbx_value(datablock_cifvalue(datablock,
                //                                         cif_tag_index, i));
                // this->data[cif_tag].push_back(pdbx_value);
            }
        }
    }
}

PDBx::~PDBx() {}

void PDBx::values(std::string cif_tag) {}
