#include "PDBx.h"

PDBx::PDBx(CIF* cif, std::vector<std::string> select_tags={}) {
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
        int64_t cif_tag_counter = 0;
        for (const std::string &cif_tag : cif_tags) {
            const ssize_t cif_tag_index =
                datablock_tag_index(datablock,
                                    const_cast<char*>(cif_tag.c_str()));

            if (cif_tag_index < 0) {
                continue;
            }

            this->cif_tag_order[cif_tag] = cif_tag_counter;
            this->in_loop[cif_tag] = are_in_loops[cif_tag_index];

            for (ssize_t i = 0; i < cif_value_lengths[cif_tag_index]; i++) {
                PDBXVALUE pdbx_value(datablock_cifvalue(datablock,
                                                        cif_tag_index, i));
                this->data[cif_tag].push_back(pdbx_value);
            }

            cif_tag_counter++;
        }
    }
}

std::vector<PDBXVALUE> PDBx::values(std::string cif_tag) {
    return this->data[cif_tag];
}

PDBXVALUE PDBx::value(std::string cif_tag, size_t index) {
    if (this->data[cif_tag].size() > 0) {
        return this->data[cif_tag][index];
    } else {
        return PDBXVALUE(".");
    }
}
