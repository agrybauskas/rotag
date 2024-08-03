#ifndef SRC_LIB_ATOMSITE_H_
#define SRC_LIB_ATOMSITE_H_

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

extern "C" {
    #include "cif.h"
    #include "cif_compiler.h"
}

#include "PDBx.h"

struct Selector {
    std::map<std::string, std::map<std::string, bool>> selection;

    Selector() {}

    Selector(std::map<std::string, std::vector<std::string>> tag_selection) {
        for (std::map<std::string, std::vector<std::string>>::iterator tag_it =
                 tag_selection.begin();
             tag_it != tag_selection.end();
             ++tag_it) {
            std::string tag_name = tag_it->first;
            std::vector<std::string> tag_values = tag_selection[tag_name];
            for (std::string& tag_value : tag_values) {
                this->selection[tag_name][tag_value] = true;
            }
        }
    }

    void add(std::string tag_name, std::string tag_value) {
        selection[tag_name][tag_value] = true;
    }
};

class AtomSite {
 private:
    enum M_TAG_INDEX {
        GROUP_PDB,              // "ATOM" or "HETATM".
        ID,                     // Atom id.
        TYPE_SYMBOL,            // Chemical element.
        LABEL_ATOM_ID,          // Atom label.
        LABEL_ALT_ID,           // Related to alt. atom position.
        LABEL_COMP_ID,          // Residue name.
        LABEL_ASYM_ID,          // Chain name.
        LABEL_ENTITY_ID,        // Molecular entity.
        LABEL_SEQ_ID,           // Residue id.
        CARTN_X,                // Cartesian x coordinates of the atom.
        CARTN_Y,                // Cartesian y coordinates of the atom.
        CARTN_Z,                // Cartesian z coordinates of the atom.
        OCCUPANCY,              // The fraction present in the site.
        B_ISO_OR_EQUIV,         // Isotropic displacement.
        AUTH_SEQ_ID,            // Author's residue id.
        AUTH_COMP_ID,           // Author's residue name.
        AUTH_ASYM_ID,           // Author's chain name.
        AUTH_ATOM_ID,           // Author's atom label.
        PDBX_PDB_MODEL_NUM,     // Model id.
        ROTAG_SELECTION_STATE,  // Marks selection state: T, S ir H.
        ROTAG_SELECTION_GROUP   // Selection group id.
    };

    const std::vector<std::string> M_TAGS = {
        "_atom_site.group_pdb",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.cartn_x",
        "_atom_site.cartn_y",
        "_atom_site.cartn_z",
        "_atom_site.occupancy",
        "_atom_site.b_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_pdb_model_num",
        "_atom_site.rotag_selection_state",
        "_atom_site.rotag_selection_group"
    };

    const M_TAG_INDEX M_ID = ID;  // Declares UNIQUE ID for the class object.

    PDBx m_data;

 public:
    AtomSite();
    explicit AtomSite(char*, bool);

    const std::vector<std::string> names();
    const std::string name(int64_t);
    std::vector<PDBXVALUE> values(std::string);
    // PDBXVALUE value(int64_t, std::string);
    // PDBXVALUE value(int64_t, int64_t);
    // void add_atom(int64_t, Atom);
    // std::vector<PDBXVALUE> ids();
    // void mark_selection(AtomSite&,
    //                     std::vector<int64_t> = {},
    //                     std::vector<int64_t> = {});
};

// AtomSite filter(AtomSite&, Selector = {{}}, Selector = {{}});

#endif  // SRC_LIB_ATOMSITE_H_
