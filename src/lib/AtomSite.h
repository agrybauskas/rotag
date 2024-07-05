#ifndef SRC_LIB_ATOMSITE_H_
#define SRC_LIB_ATOMSITE_H_

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gemmi/cif.hpp>

extern "C" {
    #include "cif.h"
    #include "cif_compiler.h"
}

#include "PDBx.h"

typedef std::map<std::string, PDBXVALUE> m_Atom;

enum M_FORMAT {
    PDBX,
    PDB
};

class AtomSite {
 private:
    typedef std::vector<PDBXVALUE> m_PDBXVALUES;
    typedef std::map<std::string, std::map<std::string, bool>> m_Selector;

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
        "group_pdb",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "cartn_x",
        "cartn_y",
        "cartn_z",
        "occupancy",
        "b_iso_or_equiv",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_pdb_model_num",
        "rotag_selection_state",
        "rotag_selection_group"
    };

    const M_TAG_INDEX M_ID = ID;  // Declares UNIQUE ID for the class object.

    std::map<int64_t, m_Atom> m_atoms = {};

 public:
    explicit AtomSite(char*, M_FORMAT = PDBX);

    const std::vector<std::string> names();
    const std::string name(int64_t);
    void mark_selection(AtomSite&, std::vector<int64_t>, std::vector<int64_t>);
    AtomSite filter(AtomSite, m_Selector, m_Selector);
    std::map<int64_t, m_Atom> atoms();
    m_Atom atom(int64_t);
    PDBXVALUE value(int64_t, int64_t);
};

#endif  // SRC_LIB_ATOMSITE_H_
