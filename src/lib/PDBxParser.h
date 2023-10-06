#ifndef SRC_LIB_PDBXPARSER_H_
#define SRC_LIB_PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>

struct Atom {
  // "_atom_site" category-related.
  std::string group_pdb;                // "ATOM" or "HETATM".
  unsigned long int id;                 // Atom id.
  std::string type_symbol;              // Chemical element.
  std::string label_atom_id;            // Atom label.
  std::string label_alt_id;             // Related to alternative atom position.
  std::string label_comp_id;            // Residue name.
  std::string label_asym_id;            // Chain name.
  unsigned long int label_entity_id;    // Molecular entity.
  unsigned long int label_seq_id;       // Residue id.
  double cartn_x, cartn_y, cartn_z;     // Cartesian coordinates of the atom.
  double occupancy;                     // The fraction present in the site.
  double b_iso_or_equiv;                // Isotropic displacement.
  unsigned long int auth_seq_id;        // Author's residue id.
  std::string auth_comp_id;             // Author's residue name.
  std::string auth_asym_id;             // Author's chain name.
  std::string auth_atom_id;             // Author's atom label.
  unsigned long int pdbx_pdb_model_num; // Model id.

  // Selection.
  std::string selection_state;          // Marks selection state: T, S ir H.
  std::string selection_group;          // Selection group id.
};

struct Selector {
  std::map<std::string, bool> group_pdb;
  std::map<unsigned long int, bool> id;
  std::map<std::string, bool> type_symbol;
  std::map<std::string, bool> label_atom_id;
  std::map<std::string, bool> label_alt_id;
  std::map<std::string, bool> label_comp_id;
  std::map<std::string, bool> label_asym_id;
  std::map<unsigned long int, bool> label_entity_id;
  std::map<unsigned long int, bool> label_seq_id;
  std::map<unsigned long int, bool> auth_seq_id;
  std::map<std::string, bool> auth_comp_id;
  std::map<std::string, bool> auth_asym_id;
  std::map<std::string, bool> auth_atom_id;
  std::map<unsigned long int, bool> pdbx_pdb_model_num;
  std::map<std::string, bool> selection_state;
  std::map<std::string, bool> selection_group;
};

typedef std::map<unsigned int, Atom> AtomSite;

AtomSite mmcif_to_atom_site(char* mmcif_file_path);

std::vector<unsigned long int>
  filter(AtomSite atom_site, Selector include={}, Selector exclude={});

AtomSite extract(AtomSite atom_site);

void mark_selection(AtomSite atom_site);

#endif  // SRC_LIB_PDBXPARSER_H_
