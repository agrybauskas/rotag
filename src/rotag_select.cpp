#include <getopt.h>
#include <iostream>
#include <string>

#include "lib/Version.h"

int main(int argc, char *argv[]) {
  const struct option longopts[] = {
    {"target",         0, 0, 't'},
    {"select",         0, 0, 's'},
    {"tags",           0, 0     },
    {"related-data",   0, 0, 'r'},
    {"pdb",            0, 0, 'p'},
    {"keep-ignored",   0, 0, 'k'},
    {"random-seed",    0, 0, 'x'},
    {"help",           0, 0, 'h'},
    {"version",        0, 0, 'v'},
    {0,                0, 0, 0  },
  };

  int index;
  int iarg = 0;

  while(iarg != -1) {
    iarg = getopt_long(argc, argv, "tsrpkxh:v:", longopts, &index);
    switch(iarg) {
      case 'h':
        std::cout << "rotag_select [options] [--] <cif-file>...\n"
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
          "        command describing the residues that will be marked\n"
          "        as target (T) (default: all).\n"
          "\n"
          "        Selection keywords (equivalent PDBx category data\n"
          "        items in parentheses):\n"
          "            all       - all atoms;\n"
          "            atomid    - atom id number;\n"
          "            atomname  - atom name;\n"
          "            atomtype  - atom symbol;\n"
          "            resid     - residue id number;\n"
          "            authresid - author residue id number;\n"
          "            resname   - residue name;\n"
          "            chain     - chain name;\n"
          "            authchain - author chain name;\n"
          "            altid     - alt atom id;\n"
          "            model     - pdbx model num;\n"
          "            mainchain - mainchain atoms;\n"
          "            sidechain - sidechain atoms;\n"
          "            hetatoms  - heteroatoms;\n"
          "            target    - target atoms (--select uses\n"
          "                        --target atoms and --target uses\n"
          "                        target atoms described in a file).\n"
          "\n"
          "        Selection operators:\n"
          "            around - atoms around <int> angstroms;\n"
          "            byres  - expands to residues;\n"
          "            rand   - selects <int> number of atoms randomly;\n"
          "            angles - assigns min <float> and max <float>\n"
          "                     angle range where selected angles should\n"
          "                     be in (e.g. angles 'chi1=90..160' or\n"
          "                     'psi=60..70, phi=50..80');\n"
          "            expand - expands the selection by one step that\n"
          "                     follows the connections.\n"
          "\n"
          "        List operators:\n"
          "            .. - range of integers (e.g. 2..5 is 2,3,4,5);\n"
          "            ,  - list of integers or keywords (e.g. 2,3 or\n"
          "                 A,B).\n"
          "\n"
          "        Map operators:\n"
          "            = - assigns value to the string (e.g. chi1=1.0).\n"
          "\n"
          "        Logical operators:\n"
          "            && - and operator;\n"
          "            || - or operator;\n"
          "            !  - negation operator;\n"
          "            () - parentheses;\n"
          "            ;  - group/order selection separator.\n"
          "\n"
          "        Annotation operator:\n"
          "            : - assigns <int> or <str> as group id.\n"
          "    -s, --select <selection-query>\n"
          "        command (same as --target) describing the atoms that\n"
          "        will be marked as selected (S) (default: target).\n"
          "    --tags <tag>[,<tag>...]\n"
          "        select PDBx tag that will be in the output\n"
          "        (default: '_atom_site,_[local]_rotamer_angle,\n"
          "         _[local]_dihedral_angle,_[local]_rotamer_energy,\n"
          "         _[local]_pairwise_energy,_[local]_energy,_[local]_rmsd').\n"
          "    -r, --related-data\n"
          "        only related data records from other categories are\n"
          "        shown when '_atom_site' records are selected.\n"
          "    -p, --pdb\n"
          "        input file is in PDB format.\n"
          "    -k, --keep-ignored\n"
          "        keep ignored atoms.\n"
          "    -x, --rand-seed <int>\n"
          "        set a seed for random (rand) selection.\n"
          "    -v, --version\n"
          "        print version" << std::endl;
        break;
      case 'v':
        std::cout << "???" << std::endl;
        break;
    }
  }

  return 0;
}
