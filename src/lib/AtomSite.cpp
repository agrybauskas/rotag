#include "AtomSite.h"

AtomSite::AtomSite() {}

AtomSite::AtomSite(char* file_path, bool is_pdb = false) {
    if (is_pdb) {
        /* gemmi::Structure structure = gemmi::read_pdb_file(file_path); */
    } else {
        cif_option_t compiler_options = cif_option_default();
        cexception_t inner;
        CIF* cif = new_cif_from_cif_file(file_path, compiler_options, &inner);

        PDBx pdbx(cif, this->M_TAGS);
        this->m_data = pdbx;

        // NOTE: needs to be refactored and ID standardised.
        for (size_t i = 0; i < this->m_data.length(this->name(ID)); i++) {
            m_id_to_index[this->m_data.value(this->name(ID), i)] = i;
        }

        delete_cif(cif);
    }
}

const std::vector<std::string> AtomSite::names() {
    return this->M_TAGS;
}

const std::string AtomSite::name(int64_t index) {
    return this->M_TAGS[index];
}

std::vector<PDBXVALUE> AtomSite::values(std::string cif_tag) {
    return this->m_data.values(cif_tag);
}

PDBXVALUE AtomSite::value(int64_t id, std::string cif_tag) {
    return this->m_data.values(cif_tag).at(this->m_id_to_index[id]);
}

PDBXVALUE AtomSite::value(int64_t id, int64_t index) {
    return this->m_data.values(this->name(index)).at(this->m_id_to_index[id]);
}

// void AtomSite::add_atom(int64_t id, Atom atom) {
//     this->m_atoms.insert(std::make_pair(id, atom));
// }

std::vector<PDBXVALUE> AtomSite::ids() {
    return this->m_data.values(this->name(M_ID));
}

// void AtomSite::mark_selection(AtomSite& atom_site,
//                               std::vector<int64_t> target_atom_ids,
//                               std::vector<int64_t> selected_atom_ids) {
//     std::map<int64_t, Atom> atoms = atom_site.atoms();
//     std::map<int64_t, Atom>::iterator atom_it;
//     for (atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
//         int64_t id = atom_it->first;
//         atom_site.m_atoms.at(id).erase("_atom_site.rotag_selection_state");
//         atom_site.m_atoms.at(id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("I"));
//     }
//     for (const int64_t &selected_atom_id : selected_atom_ids) {
//         atom_site.m_atoms.at(selected_atom_id).erase(
//             "_atom_site.rotag_selection_state");
//         atom_site.m_atoms.at(selected_atom_id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("S"));
//     }
//     for (const int64_t &target_atom_id : target_atom_ids) {
//         atom_site.m_atoms.at(target_atom_id).erase(
//             "_atom_site.rotag_selection_state");
//         atom_site.m_atoms.at(target_atom_id).emplace(
//             "_atom_site.rotag_selection_state", PDBXVALUE("T"));
//     }
// }

AtomSite filter(AtomSite& atom_site,
                Selector include,
                Selector exclude) {
    if (atom_site.ids().empty()) {
        /* TODO(algirdas): Error or warning. Message: no atoms were loaded from
           "_atom_site".*/
    }

    AtomSite filtered_atom_site;
//     std::map<int64_t, Atom> atoms = atom_site.atoms();
//     std::map<int64_t, Atom>::iterator atom_it;
//     for (atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
//         int64_t id = atom_it->first;
//         Atom atom = atom_site.atom(id);
//         bool keep_atom = false;
//         for (Atom::iterator tag_it = atom.begin();
//              tag_it != atom.end();
//              ++tag_it) {
//             std::string cif_tag = tag_it->first;
//             std::string value = atom_site.value(id, cif_tag);

//             if (!include.selection[cif_tag].empty() &&
//                 include.selection[cif_tag][value]) {
//                 keep_atom = true;
//                 break;
//             }

//             if (!exclude.selection[cif_tag].empty() &&
//                 exclude.selection[cif_tag][value]) {
//                 keep_atom = false;
//                 break;
//             }
//         }

//         if (!keep_atom) {
//             continue;
//         }

//         filtered_atom_site.add_atom(id, atom_site.atom(id));
//     }

    return filtered_atom_site;
}
