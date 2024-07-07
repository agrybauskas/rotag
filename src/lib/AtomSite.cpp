#include "AtomSite.h"

AtomSite::AtomSite() {}

AtomSite::AtomSite(char* file_path, bool is_pdb = false) {
    if (is_pdb) {
        gemmi::Structure structure = gemmi::read_pdb_file(file_path);
    } else {
        cif_option_t compiler_options = cif_option_default();
        cexception_t inner;
        CIF* cif = new_cif_from_cif_file(file_path, compiler_options, &inner);

        PDBx pdbx(cif, this->M_TAGS);

        delete_cif(cif);

        std::vector<PDBXVALUE> id_values = pdbx.values(this->name(M_ID));
        for (int m_tag_index = GROUP_PDB;
             m_tag_index <= ROTAG_SELECTION_GROUP;
             m_tag_index++) {
            for (size_t i = 0; i < id_values.size(); i++) {
                this->m_atoms[id_values[i]].insert(
                    std::make_pair(this->name(m_tag_index),
                                   pdbx.value(this->name(m_tag_index), i)));
            }
        }
    }
}

const std::vector<std::string> AtomSite::names() {
    return this->M_TAGS;
}

const std::string AtomSite::name(int64_t index) {
    return this->M_TAGS[index];
}

std::map<int64_t, Atom> AtomSite::atoms() {
    return this->m_atoms;
}

Atom AtomSite::atom(int64_t id) {
    return this->m_atoms.at(id);
}

std::vector<PDBXVALUE> AtomSite::values(std::string cif_tag) {
    std::vector<PDBXVALUE> values;
    std::map<int64_t, Atom> atoms = this->atoms();
    std::map<int64_t, Atom>::iterator atom_it;
    for (atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
        int64_t id = atom_it->first;
        values.push_back(this->value(id, cif_tag));
    }
    return values;
}

PDBXVALUE AtomSite::value(int64_t id, std::string cif_tag) {
    return this->atom(id).at(cif_tag);
}

PDBXVALUE AtomSite::value(int64_t id, int64_t index) {
    return this->atom(id).at(this->name(index));
}

void AtomSite::add_atom(int64_t id, Atom atom) {
    this->m_atoms.insert(std::make_pair(id, atom));
}

std::vector<PDBXVALUE> AtomSite::ids() {
    return this->values(this->name(M_ID));
}

void AtomSite::mark_selection(AtomSite& atom_site,
                              std::vector<int64_t> target_atom_ids,
                              std::vector<int64_t> selected_atom_ids) {
    std::map<int64_t, Atom> atoms = atom_site.atoms();
    std::map<int64_t, Atom>::iterator atom_it;
    for (atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
        int64_t id = atom_it->first;
        atom_site.m_atoms.at(id).erase("_atom_site.rotag_selection_state");
        atom_site.m_atoms.at(id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("I"));
    }
    for (const int64_t &selected_atom_id : selected_atom_ids) {
        atom_site.m_atoms.at(selected_atom_id).erase(
            "_atom_site.rotag_selection_state");
        atom_site.m_atoms.at(selected_atom_id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("S"));
    }
    for (const int64_t &target_atom_id : target_atom_ids) {
        atom_site.m_atoms.at(target_atom_id).erase(
            "_atom_site.rotag_selection_state");
        atom_site.m_atoms.at(target_atom_id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("T"));
    }
}

AtomSite filter(AtomSite& atom_site,
                Selector include,
                Selector exclude) {
    if (atom_site.atoms().empty()) {
        /* TODO(algirdas): Error or warning. Message: no atoms were loaded from
           "_atom_site".*/
    }

    AtomSite filtered_atom_site;
    std::map<int64_t, Atom> atoms = atom_site.atoms();
    std::map<int64_t, Atom>::iterator atom_it;
    for (atom_it = atoms.begin(); atom_it != atoms.end(); ++atom_it) {
        int64_t id = atom_it->first;
        Atom atom = atom_site.atom(id);
        bool keep_atom = true;
        for (Atom::iterator tag_it = atom.begin();
             tag_it != atom.end();
             ++tag_it) {
            std::string cif_tag = tag_it->first;
            std::string value = atom_site.value(id, cif_tag);

            if (!include[cif_tag].empty() && !include[cif_tag][value]) {
                keep_atom = false;
                break;
            }

            if (!exclude[cif_tag].empty() && exclude[cif_tag][value]) {
                keep_atom = false;
                break;
            }
        }

        if (!keep_atom) {
            continue;
        }

        filtered_atom_site.add_atom(id, atom_site.atom(id));
    }

    return filtered_atom_site;
}
