#ifndef SRC_LIB_PDBXPARSER_H_
#define SRC_LIB_PDBXPARSER_H_

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

class AtomSite {
 private:

 public:
    AtomSite(char*) {};
};

// enum ATOM_SITE_INDEX {
//     GROUP_PDB,             // "ATOM" or "HETATM".
//     ID,                    // Atom id.
//     TYPE_SYMBOL,           // Chemical element.
//     LABEL_ATOM_ID,         // Atom label.
//     LABEL_ALT_ID,          // Related to alt. atom position.
//     LABEL_COMP_ID,         // Residue name.
//     LABEL_ASYM_ID,         // Chain name.
//     LABEL_ENTITY_ID,       // Molecular entity.
//     LABEL_SEQ_ID,          // Residue id.
//     CARTN_X,               // Cartesian x coordinates of the atom.
//     CARTN_Y,               // Cartesian y coordinates of the atom.
//     CARTN_Z,               // Cartesian z coordinates of the atom.
//     OCCUPANCY,             // The fraction present in the site.
//     B_ISO_OR_EQUIV,        // Isotropic displacement.
//     AUTH_SEQ_ID,           // Author's residue id.
//     AUTH_COMP_ID,          // Author's residue name.
//     AUTH_ASYM_ID,          // Author's chain name.
//     AUTH_ATOM_ID,          // Author's atom label.
//     PDBX_PDB_MODEL_NUM,    // Model id.
//     ROTAG_SELECTION_STATE, // Marks selection state: T, S ir H.
//     ROTAG_SELECTION_GROUP  // Selection group id.
// };

// struct ATOM_SITE_TAGS {
//     const std::vector<std::string> NAMES = {
//         "_atom_site.group_pdb",
//         "_atom_site.id",
//         "_atom_site.type_symbol",
//         "_atom_site.label_atom_id",
//         "_atom_site.label_alt_id",
//         "_atom_site.label_comp_id",
//         "_atom_site.label_asym_id",
//         "_atom_site.label_entity_id",
//         "_atom_site.label_seq_id",
//         "_atom_site.cartn_x",
//         "_atom_site.cartn_y",
//         "_atom_site.cartn_z",
//         "_atom_site.occupancy",
//         "_atom_site.b_iso_or_equiv",
//         "_atom_site.auth_seq_id",
//         "_atom_site.auth_comp_id",
//         "_atom_site.auth_asym_id",
//         "_atom_site.auth_atom_id",
//         "_atom_site.pdbx_pdb_model_num",
//         "_atom_site.rotag_selection_state",
//         "_atom_site.rotag_selection_group"
//     };

//     const std::string name(int index) {
//         return NAMES[index];
//     }

//     const std::vector<std::string> names() {
//         return NAMES;
//     }
// };

// typedef std::map<std::string, PDBXVALUE> Atom;
// typedef std::vector<PDBXVALUE> PDBXVALUES;
// typedef std::map<int64_t, Atom> AtomSite;
// typedef std::map<std::string, std::map<std::string, bool>> Selector;

// AtomSite pdbx_to_atom_site(char*);

// AtomSite filter(AtomSite, Selector, Selector);

// void mark_selection(AtomSite&, std::vector<int64_t>, std::vector<int64_t>);

#endif  // SRC_LIB_PDBXPARSER_H_
