#include "PDBxParser.h"

AtomSite pdbx_to_atom_site(char* pdbx_file_path) {
    AtomSite atom_site = {};

    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    CIF* cif = new_cif_from_cif_file(pdbx_file_path, compiler_options, &inner);

    PDBx pdbx(cif, ATOM_SITE_TAGS);

    delete_cif(cif);

    for (size_t i = 0; i < pdbx.values("_atom_site.id").size(); i++) {
        int64_t id = pdbx.values("_atom_site.id")[i];
        Atom atom = {};
        for (const std::string &cif_tag : ATOM_SITE_TAGS) {
            if (pdbx.values(cif_tag).size() > 0) {
                atom.insert(std::make_pair(cif_tag, pdbx.values(cif_tag)[i]));
            }
        }
        atom_site.insert(std::make_pair(id, atom));
    }

    return atom_site;
}

AtomSite filter(AtomSite atom_site, Selector include, Selector exclude) {
    if (atom_site.empty()) {
        /* TODO(algirdas): Error or warning. Message: no atoms were loaded from
           "_atom_site".*/
    }

    AtomSite filtered_atom_site = {};
    for (AtomSite::iterator it_i = atom_site.begin(); it_i != atom_site.end(); ++it_i) {
        int64_t id = it_i->first;
        Atom atom = atom_site.at(id);
        bool keep_atom = true;
        for (Atom::iterator it_j = atom.begin(); it_j != atom.end(); ++it_j) {
            std::string cif_tag = it_j->first;
            std::string value = atom.at(cif_tag);

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

        filtered_atom_site[id] = atom_site[id];
    }

    return filtered_atom_site;
}

void mark_selection(AtomSite atom_site,
                    std::vector<int64_t> target_atom_ids,
                    std::vector<int64_t> selected_atom_ids) {
    for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
        int64_t id = it->first;
        std::make_pair(atom_site[id], PDBXVALUE("_atom_site.rotag_selection_state", "I"));
    }
    for (const int64_t &selected_atom_id : selected_atom_ids) {
        /* atom_site[selected_atom_id] = PDBXVALUE("_atom_site.rotag_selection_state", "S"); */
    }
    for (const int64_t &target_atom_id : target_atom_ids) {
        /* atom_site[target_atom_id] = PDBXVALUE("_atom_site.rotag_selection_state", "T"); */
    }
}
