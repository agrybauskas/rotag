#include "PDBxParser.h"

// AtomSite pdbx_to_atom_site(char* pdbx_file_path) {
//     AtomSite atom_site = {};

//     cif_option_t compiler_options = cif_option_default();
//     cexception_t inner;
//     CIF* cif = new_cif_from_cif_file(pdbx_file_path, compiler_options, &inner);

//     PDBx pdbx(cif, ATOM_SITE_TAGS().names());

//     delete_cif(cif);

//     ATOM_SITE_TAGS TAGS;
//     PDBXVALUES group_pdbs = pdbx.values(TAGS.name(GROUP_PDB));
//     PDBXVALUES atom_ids = pdbx.values(TAGS.name(ID));
//     PDBXVALUES type_symbols = pdbx.values(TAGS.name(TYPE_SYMBOL));
//     PDBXVALUES label_atom_ids = pdbx.values(TAGS.name(LABEL_ATOM_ID));
//     PDBXVALUES label_alt_ids = pdbx.values(TAGS.name(LABEL_ALT_ID));
//     PDBXVALUES label_comp_ids = pdbx.values(TAGS.name(LABEL_COMP_ID));
//     PDBXVALUES label_asym_ids = pdbx.values(TAGS.name(LABEL_ASYM_ID));
//     PDBXVALUES label_entity_ids = pdbx.values(TAGS.name(LABEL_ENTITY_ID));
//     PDBXVALUES label_seq_ids = pdbx.values(TAGS.name(LABEL_SEQ_ID));
//     PDBXVALUES cartn_xs = pdbx.values(TAGS.name(CARTN_X));
//     PDBXVALUES cartn_ys = pdbx.values(TAGS.name(CARTN_Y));
//     PDBXVALUES cartn_zs = pdbx.values(TAGS.name(CARTN_Z));
//     PDBXVALUES occupancies = pdbx.values(TAGS.name(OCCUPANCY));
//     PDBXVALUES b_iso_or_equivs = pdbx.values(TAGS.name(B_ISO_OR_EQUIV));
//     PDBXVALUES auth_seq_ids = pdbx.values(TAGS.name(AUTH_SEQ_ID));
//     PDBXVALUES auth_comp_ids = pdbx.values(TAGS.name(AUTH_COMP_ID));
//     PDBXVALUES auth_asym_ids = pdbx.values(TAGS.name(AUTH_ASYM_ID));
//     PDBXVALUES auth_atom_ids = pdbx.values(TAGS.name(AUTH_ATOM_ID));
//     PDBXVALUES pdbx_pdb_model_nums = pdbx.values(TAGS.name(PDBX_PDB_MODEL_NUM));
//     PDBXVALUES rotag_selection_states =
//         pdbx.values(TAGS.name(ROTAG_SELECTION_STATE));
//     PDBXVALUES rotag_selection_groups =
//         pdbx.values(TAGS.name(ROTAG_SELECTION_GROUP));

//     for (size_t i = 0; i < atom_ids.size(); i++) {
//         Atom atom = {
//             {TAGS.name(GROUP_PDB), group_pdbs[i]},
//             {TAGS.name(ID), atom_ids[i]},
//             {TAGS.name(TYPE_SYMBOL), type_symbols[i]},
//             {TAGS.name(LABEL_ATOM_ID), label_atom_ids[i]},
//             {TAGS.name(LABEL_ALT_ID), label_alt_ids[i]},
//             {TAGS.name(LABEL_COMP_ID), label_comp_ids[i]},
//             {TAGS.name(LABEL_ASYM_ID), label_asym_ids[i]},
//             {TAGS.name(LABEL_ENTITY_ID), label_entity_ids[i]},
//             {TAGS.name(LABEL_SEQ_ID), label_seq_ids[i]},
//             {TAGS.name(CARTN_X), cartn_xs[i]},
//             {TAGS.name(CARTN_Y), cartn_ys[i]},
//             {TAGS.name(CARTN_Z), cartn_zs[i]},
//             {TAGS.name(OCCUPANCY), occupancies[i]},
//             {TAGS.name(B_ISO_OR_EQUIV), b_iso_or_equivs[i]},
//             {TAGS.name(AUTH_SEQ_ID), auth_seq_ids[i]},
//             {TAGS.name(AUTH_COMP_ID), auth_comp_ids[i]},
//             {TAGS.name(AUTH_ASYM_ID), auth_asym_ids[i]},
//             {TAGS.name(AUTH_ATOM_ID), auth_atom_ids[i]},
//             {TAGS.name(PDBX_PDB_MODEL_NUM), pdbx_pdb_model_nums[i]},
//             // {TAGS.name(ROTAG_SELECTION_STATE), rotag_selection_states[i]},
//             // {TAGS.name(ROTAG_SELECTION_GROUP), rotag_selection_groups[i]}
//         };
//         atom_site.insert(std::make_pair(atom_ids[i], atom));
//     }

//     return atom_site;
// }

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
