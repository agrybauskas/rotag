#include "AtomSite.h"

AtomSite::AtomSite(char* pdbx_file_path) {
    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    CIF* cif = new_cif_from_cif_file(pdbx_file_path, compiler_options, &inner);

    PDBx pdbx(cif, this->M_TAGS);

    delete_cif(cif);

    m_PDBXVALUES id_values = pdbx.values(this->name(ID));
    for (int m_tag_index = GROUP_PDB;
         m_tag_index <= ROTAG_SELECTION_GROUP;
         m_tag_index++) {
        for (size_t i = 0; i < id_values.size(); i++) {
            this->m_atoms[id_values[i]].insert(
                std::make_pair(this->name(m_tag_index),
                               pdbx.value(this->M_TAGS[m_tag_index], i)));
        }
    }
}

const std::vector<std::string> AtomSite::names() {
    return this->M_TAGS;
}

const std::string AtomSite::name(int index) {
    return this->M_TAGS[index];
}

std::map<int64_t, m_Atom> AtomSite::atoms() {
    return this->m_atoms;
}

m_Atom AtomSite::atom(int64_t id) {
    return this->atom(id);
}

// AtomSite filter(AtomSite atom_site,
//                 Selector include = {{}},
//                 Selector exclude = {{}}) {
//     if (atom_site.empty()) {
//         /* TODO(algirdas): Error or warning. Message: no atoms were loaded from
//            "_atom_site".*/
//     }

//     AtomSite filtered_atom_site = {};
//     for (AtomSite::iterator it_i = atom_site.begin(); it_i != atom_site.end(); ++it_i) {
//         int64_t id = it_i->first;
//         Atom atom = atom_site.at(id);
//         bool keep_atom = true;
//         for (Atom::iterator it_j = atom.begin(); it_j != atom.end(); ++it_j) {
//             std::string cif_tag = it_j->first;
//             std::string value = atom.at(cif_tag);

//             if (!include[cif_tag].empty() && !include[cif_tag][value]) {
//                 keep_atom = false;
//                 break;
//             }

//             if (!exclude[cif_tag].empty() && exclude[cif_tag][value]) {
//                 keep_atom = false;
//                 break;
//             }
//         }

//         if (!keep_atom) {
//             continue;
//         }

//         filtered_atom_site[id] = atom_site[id];
//     }

//     return filtered_atom_site;
// }

// void mark_selection(AtomSite& atom_site,
//                     std::vector<int64_t> target_atom_ids = {},
//                     std::vector<int64_t> selected_atom_ids = {}) {
//     for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
//         int64_t id = it->first;
//         atom_site.at(id).erase("_atom_site.rotag_selection_state");
//         atom_site.at(id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("I"));
//     }
//     for (const int64_t &selected_atom_id : selected_atom_ids) {
//         atom_site.at(selected_atom_id).erase("_atom_site.rotag_selection_state");
//         atom_site.at(selected_atom_id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("S"));
//     }
//     for (const int64_t &target_atom_id : target_atom_ids) {
//         atom_site.at(target_atom_id).erase("_atom_site.rotag_selection_state");
//         atom_site.at(target_atom_id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("T"));
//     }
// }
