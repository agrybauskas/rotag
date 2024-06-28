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

const std::vector<std::string> ATOM_SITE_TAGS = {
    // "_atom_site" category-related.
    "_atom_site.group_pdb",              // "ATOM" or "HETATM".
    "_atom_site.id",                     // Atom id.
    "_atom_site.type_symbol",            // Chemical element.
    "_atom_site.label_atom_id",          // Atom label.
    "_atom_site.label_alt_id",           // Related to alt. atom position.
    "_atom_site.label_comp_id",          // Residue name.
    "_atom_site.label_asym_id",          // Chain name.
    "_atom_site.label_entity_id",        // Molecular entity.
    "_atom_site.label_seq_id",           // Residue id.
    "_atom_site.cartn_x",                // Cartesian x coordinates of the atom.
    "_atom_site.cartn_y",                // Cartesian y coordinates of the atom.
    "_atom_site.cartn_z",                // Cartesian z coordinates of the atom.
    "_atom_site.occupancy",              // The fraction present in the site.
    "_atom_site.b_iso_or_equiv",         // Isotropic displacement.
    "_atom_site.auth_seq_id",            // Author's residue id.
    "_atom_site.auth_comp_id",           // Author's residue name.
    "_atom_site.auth_asym_id",           // Author's chain name.
    "_atom_site.auth_atom_id",           // Author's atom label.
    "_atom_site.pdbx_pdb_model_num",     // Model id.

    // Selection that is specific to rotag.
    "_atom_site.rotag_selection_state",  // Marks selection state: T, S ir H.
    "_atom_site.rotag_selection_group"   // Selection group id.
};

typedef std::map<std::string, PDBXVALUE> Atom;
typedef std::map<int64_t, Atom> AtomSite;
typedef std::map<std::string, std::map<std::string, bool>> Selector;

AtomSite pdbx_to_atom_site(char* pdbx_file_path);

AtomSite filter(AtomSite atom_site,
                Selector include = {{}},
                Selector exclude = {{}});

void mark_selection(AtomSite& atom_site,
                    std::vector<int64_t> target_atom_ids={},
                    std::vector<int64_t> selected_atom_ids={});

#endif  // SRC_LIB_PDBXPARSER_H_
