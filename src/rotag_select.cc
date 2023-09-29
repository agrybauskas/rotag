#include <getopt.h>

#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>

#include "lib/ForceField/Parameters.h"
#include "lib/PDBxParser.h"
#include "lib/Version.h"

extern char *optarg;
extern int optind, opterr, optopt;

char* progname;

int main(int argc, char *argv[]) {
  progname = argv[0];

  // Defaults.
  std::string target_cmds = "all";
  std::string select_cmds = "target";
  std::string tags =
    "_atom_site,_[local]_rotamer_angle,_[local]_dihedral_angle,"
    "_[local]_rotamer_energy,_[local]_pairwise_energy,_[local]_energy,"
    "_[local]_rmsd";
  bool is_related = false;
  bool is_pdb = false;
  bool keep_ignored = false;
  int random_seed = 23;
  std::vector<std::string> category_list = {
    "_atom_site", "_[local]_rotamer_angle", "_[local]_dihedral_angle",
    "_[local]_rotamer_energy", "_[local]_pairwise_energy", "_[local]_energy",
    "_[local]_rmsd"
  };

  const struct option longopts[] = {
    {"target",         1, 0, 't'},
    {"select",         1, 0, 's'},
    {"tags",           1, 0, 0  },
    {"related-data",   0, 0, 'r'},
    {"pdb",            0, 0, 'p'},
    {"keep-ignored",   0, 0, 'k'},
    {"random-seed",    1, 0, 'x'},
    {"help",           0, 0, 'h'},
    {"version",        0, 0, 'v'},
    {0,                0, 0, 0  },
  };

  int index;
  int iarg = 0;

  while(iarg != -1) {
    iarg = getopt_long(argc, argv, "t:s:0:rpkx:hv", longopts, &index);
    switch(iarg) {
      case 't':
        target_cmds = optarg;
        break;
      case 's':
        select_cmds = optarg;
        break;
      case 0:
        tags = optarg;
        break;
      case 'r':
        is_related = true;
        break;
      case 'p':
        is_pdb = true;
        break;
      case 'k':
        keep_ignored = true;
        break;
      case 'x':
        random_seed = atoi(optarg);
        break;
      case 'h':
        std::cout <<
"rotag_select [options] [--] <cif-file>...\n"
"    select and mark atoms of interest by adding selection state [T|S|I] to\n"
"    _atom_site category in PDBx.\n"
"\n"
"Usage:\n"
"    rotag_select \\\n"
"        --target 'resname ASP' --select 'resid 1..10' input.cif > output.cif\n"
"    rotag_select \\\n"
"        --target 'model 5' \\\n"
"        --target 'model 10' input1.cif input2.cif > output.cif\n"
"\n"
"Options:\n"
"    -t, --target <selection-query>\n"
"                        command describing the residues that will be marked\n"
"                        as target (T) (default: all).\n"
"\n"
"                        Selection keywords (equivalent PDBx category data\n"
"                        items in parentheses):\n"
"                            all       - all atoms;\n"
"                            atomid    - atom id number;\n"
"                            atomname  - atom name;\n"
"                            atomtype  - atom symbol;\n"
"                            resid     - residue id number;\n"
"                            authresid - author residue id number;\n"
"                            resname   - residue name;\n"
"                            chain     - chain name;\n"
"                            authchain - author chain name;\n"
"                            altid     - alt atom id;\n"
"                            model     - pdbx model num;\n"
"                            mainchain - mainchain atoms;\n"
"                            sidechain - sidechain atoms;\n"
"                            hetatoms  - heteroatoms;\n"
"                            target    - target atoms (--select uses\n"
"                                        --target atoms and --target uses\n"
"                                        target atoms described in a file).\n"
"\n"
"                        Selection operators:\n"
"                            around - atoms around <int> angstroms;\n"
"                            byres  - expands to residues;\n"
"                            rand   - selects <int> number of atoms randomly;\n"
"                            angles - assigns min <float> and max <float>\n"
"                                     angle range where selected angles should\n"
"                                     be in (e.g. angles 'chi1=90..160' or\n"
"                                     'psi=60..70, phi=50..80');\n"
"                            expand - expands the selection by one step that\n"
"                                     follows the connections.\n"
"\n"
"                        List operators:\n"
"                            .. - range of integers (e.g. 2..5 is 2,3,4,5);\n"
"                            ,  - list of integers or keywords (e.g. 2,3 or\n"
"                                 A,B).\n"
"\n"
"                        Map operators:\n"
"                            = - assigns value to the string (e.g. chi1=1.0).\n"
"\n"
"                        Logical operators:\n"
"                            && - and operator;\n"
"                            || - or operator;\n"
"                            !  - negation operator;\n"
"                            () - parentheses;\n"
"                            ;  - group/order selection separator.\n"
"\n"
"                        Annotation operator:\n"
"                            : - assigns <int> or <str> as group id.\n"
"    -s, --select <selection-query>\n"
"                        command (same as --target) describing the atoms that\n"
"                        will be marked as selected (S) (default: target).\n"
"    --tags <tag>[,<tag>...]\n"
"                        select PDBx tag that will be in the output\n"
"                        (default: '_atom_site,\n"
"                                   _[local]_rotamer_angle,\n"
"                                   _[local]_dihedral_angle,\n"
"                                   _[local]_rotamer_energy,\n"
"                                   _[local]_pairwise_energy,\n"
"                                   _[local]_energy,\n"
"                                   _[local]_rmsd').\n"
"    -r, --related-data\n"
"                        only related data records from other categories are\n"
"                        shown when '_atom_site' records are selected.\n"
"    -p, --pdb\n"
"                        input file is in PDB format.\n"
"    -k, --keep-ignored\n"
"                        keep ignored atoms.\n"
"    -x, --rand-seed <int>\n"
"                        set a seed for random (rand) selection.\n"
"    -v, --version\n"
"                        print version" << std::endl;
        break;
      case 'v':
        std::cout << version() << std::endl;
        break;
    }
  }

  std::vector<std::string> tags_list;
  boost::replace_all(tags, " ", "");
  boost::split(tags_list, tags, boost::is_any_of(","));

  Parameters parameters(argv[0]);

  std::cout << parameters.HYDROGEN_NAMES["ASP"]["CA"][0] << std::endl;

  return 0;
}
