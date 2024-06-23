#include "PDBxParser.h"

#include <iostream>

AtomSite pdbx_to_atom_site(char* pdbx_file_path) {
    AtomSite atom_site = {};

    cif_option_t compiler_options = cif_option_default();
    cexception_t inner;
    PDBx pdbx(new_cif_from_cif_file(pdbx_file_path, compiler_options, &inner),
              ATOM_SITE_TAGS);

    for (size_t i = 0; i < pdbx.values("_atom_site.id").size(); i++) {
        std::string id = pdbx.values("_atom_site.id")[i];
        for (const std::string &cif_tag : ATOM_SITE_TAGS) {
            // if (cif_tag_index_lookup_table[cif_tag] > 0) {
            //     atom_site[id][cif_tag] = value_scalar(datablock_cifvalue(
            //         datablock, cif_tag_index_lookup_table[cif_tag], i));
            // }
        }
    }

    return atom_site;
}

// AtomSite filter(AtomSite atom_site, Selector include, Selector exclude) {
//   if (atom_site.empty()) {
//     // TODO: Error or warning. Message: no atoms were loaded from "_atom_site".
//   }

//   AtomSite filtered_atom_site = {};
//   for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
//     std::string id = it->first;
//     bool keep_atom = true;
//     for (const std::string &cif_tag : ATOM_SITE_TAGS) {
//       std::string value = atom_site[id][cif_tag];
//       if (!include[cif_tag].empty() && !include[cif_tag][value]) {
//         keep_atom = false;
//         break;
//       }

//       if (!exclude[cif_tag].empty() && exclude[cif_tag][value]) {
//         keep_atom = false;
//         break;
//       }
//     }

//     if (!keep_atom) {
//       continue;
//     }

//     filtered_atom_site[id] = atom_site[id];
//   }

//   return filtered_atom_site;
// }

// std::vector<std::vector<std::string>>
// extract(AtomSite atom_site, std::vector<std::string> cif_tags ) {
//   std::vector<std::vector<std::string>> atoms_data = {};
//   for (AtomSite::iterator it = atom_site.begin(); it != atom_site.end(); ++it) {
//     std::string id = it->first;
//     std::vector<std::string> atom_data = {};
//     for (const std::string &cif_tag : cif_tags) {
//       atom_data.push_back(atom_site[id][cif_tag]);
//     }
//     atoms_data.push_back(atom_data);
//   }
//   return atoms_data;
// }

// void mark_selection(
//   AtomSite* atom_site,
//   std::vector<std::string> target_atom_ids,
//   std::vector<std::string> selected_atom_ids) {
//   for (AtomSite::iterator it = atom_site->begin();
//        it != atom_site->end();
//        ++it) {
//     std::string id = it->first;
//     (*atom_site)[id]["_atom_site.rotag_selection_state"] = "I";
//   }

//   for (const std::string &selected_atom_id : selected_atom_ids) {
//     (*atom_site)[selected_atom_id]["_atom_site.rotag_selection_state"] = "S";
//   }

//   for (const std::string &target_atom_id : target_atom_ids) {
//     (*atom_site)[target_atom_id]["_atom_site.rotag_selection_state"] = "T";
//   }

//   return;
// }
