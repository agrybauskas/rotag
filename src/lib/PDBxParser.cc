#include "PDBxParser.h"

AtomSite pdbx_to_atom_site(char* pdbx_file_path) {
    AtomSite atom_site = {};

    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    CIF* cif = new_cif_from_cif_file(pdbx_file_path, compiler_options, &inner);

    PDBx pdbx(cif, ATOM_SITE_TAGS);

    delete_cif(cif);

    // TODO: refactor with enums or similar.
    PDBXVALUES group_pdbs = pdbx.values("_atom_site.group_pdb");
    PDBXVALUES atom_ids = pdbx.values("_atom_site.id");
    PDBXVALUES type_symbols = pdbx.values("_atom_site.type_symbol");
    PDBXVALUES label_atom_ids = pdbx.values("_atom_site.label_atom_id");
    PDBXVALUES label_alt_ids = pdbx.values("_atom_site.label_alt_id");
    PDBXVALUES label_comp_ids = pdbx.values("_atom_site.label_comp_id");
    PDBXVALUES label_asym_ids = pdbx.values("_atom_site.label_asym_id");
    PDBXVALUES label_entity_ids = pdbx.values("_atom_site.label_entity_id");
    PDBXVALUES label_seq_ids = pdbx.values("_atom_site.label_seq_id");
    PDBXVALUES cartn_xs = pdbx.values("_atom_site.cartn_x");
    PDBXVALUES cartn_ys = pdbx.values("_atom_site.cartn_y");
    PDBXVALUES cartn_zs = pdbx.values("_atom_site.cartn_z");
    PDBXVALUES occupancies = pdbx.values("_atom_site.occupancy");
    PDBXVALUES b_iso_or_equivs = pdbx.values("_atom_site.b_iso_or_equiv");
    PDBXVALUES auth_seq_ids = pdbx.values("_atom_site.auth_seq_id");
    PDBXVALUES auth_comp_ids = pdbx.values("_atom_site.auth_comp_id");
    PDBXVALUES auth_asym_ids = pdbx.values("_atom_site.auth_asym_id");
    PDBXVALUES auth_atom_ids = pdbx.values("_atom_site.auth_atom_id");
    PDBXVALUES pdbx_pdb_model_nums =
        pdbx.values("_atom_site.pdbx_pdb_model_num");
    PDBXVALUES rotag_selection_states =
        pdbx.values("_atom_site.rotag_selection_state");
    PDBXVALUES rotag_selection_groups =
        pdbx.values("_atom_site.rotag_selection_group");

    for (size_t i = 0; i < atom_ids.size(); i++) {
        Atom atom = {
            {"_atom_site.group_pdb", group_pdbs[i]},
            {"_atom_site.id", atom_ids[i]},
            {"_atom_site.type_symbol", type_symbols[i]},
            {"_atom_site.label_atom_id", label_atom_ids[i]},
            {"_atom_site.label_alt_id", label_alt_ids[i]},
            {"_atom_site.label_comp_id", label_comp_ids[i]},
            {"_atom_site.label_asym_id", label_asym_ids[i]},
            {"_atom_site.label_entity_id", label_entity_ids[i]},
            {"_atom_site.label_seq_id", label_seq_ids[i]},
            {"_atom_site.cartn_x", cartn_xs[i]},
            {"_atom_site.cartn_y", cartn_ys[i]},
            {"_atom_site.cartn_z", cartn_zs[i]},
            {"_atom_site.occupancy", occupancies[i]},
            {"_atom_site.b_iso_or_equiv", b_iso_or_equivs[i]},
            {"_atom_site.auth_seq_id", auth_seq_ids[i]},
            {"_atom_site.auth_comp_id", auth_comp_ids[i]},
            {"_atom_site.auth_asym_id", auth_asym_ids[i]}
        };
        atom_site.insert(std::make_pair(atom_ids[i], atom));
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

void mark_selection(AtomSite& atom_site,
                    std::vector<int64_t> target_atom_ids,
                    std::vector<int64_t> selected_atom_ids) {
    for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
        int64_t id = it->first;
        atom_site.at(id).erase("_atom_site.rotag_selection_state");
        atom_site.at(id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("I"));
    }
    for (const int64_t &selected_atom_id : selected_atom_ids) {
        atom_site.at(selected_atom_id).erase("_atom_site.rotag_selection_state");
        atom_site.at(selected_atom_id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("S"));
    }
    for (const int64_t &target_atom_id : target_atom_ids) {
        atom_site.at(target_atom_id).erase("_atom_site.rotag_selection_state");
        atom_site.at(target_atom_id).emplace(
            "_atom_site.rotag_selection_state", PDBXVALUE("T"));
    }
}
