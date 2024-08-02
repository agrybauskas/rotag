#include <getopt.h>

#include <cstdio>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>

#include "lib/AtomSite.h"
#include "lib/Combinatorics.h"
#include "lib/ForceField/Parameters.h"
#include "lib/Grammar/SelectionParser.h"
#include "lib/Version.h"

extern char *optarg;
extern int optind, opterr, optopt;

char* progname;

int main(int argc, char *argv[]) {
    progname = argv[0];

    // Defaults.
    std::string target_cmd = "all";
    std::string select_cmd = "target";
    std::string tags =
        "_atom_site,_rotag_rotamer_angle,_rotag_dihedral_angle,"
        "_rotag_rotamer_energy,_rotag_pairwise_energy,_rotag_energy,"
        "_rotag_rmsd";
    bool is_related = false;
    bool is_pdb = false;
    bool keep_ignored = false;
    int64_t random_seed = 23;
    std::vector<std::string> category_list = {
        "_atom_site", "_rotag_rotamer_angle", "_rotag_dihedral_angle",
        "_rotag_rotamer_energy", "_rotag_pairwise_energy", "_rotag_energy",
        "_rotag_rmsd"
    };

    const struct option longopts[] = {
        {"target",       1, 0, 't'},
        {"select",       1, 0, 's'},
        {"tags",         1, 0, 0  },
        {"related-data", 0, 0, 'r'},
        {"pdb",          0, 0, 'p'},
        {"keep-ignored", 0, 0, 'k'},
        {"random-seed",  1, 0, 'x'},
        {"help",         0, 0, 'h'},
        {"version",      0, 0, 'v'},
        {0,              0, 0, 0  },
    };

    int index;
    int iarg = 0;

    while (iarg != -1) {
        iarg = getopt_long(argc, argv, "t:s:0:rpkx:hv", longopts, &index);
        switch (iarg) {
            case 't':
                target_cmd = optarg;
                break;
            case 's':
                select_cmd = optarg;
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
"        command describing the residues that will be marked as target (T)\n"
"        (default: all).\n"
"\n"
"        Selection keywords (equivalent PDBx category dataitems in\n"
"        parentheses):\n"
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
"            target    - target atoms (--select uses --target atoms and\n"
"                        --target uses target atoms described in a file).\n"
"\n"
"        Selection operators:\n"
"            around - atoms around <int> angstroms;\n"
"            byres  - expands to residues;\n"
"            rand   - selects <int> number of atoms randomly;\n"
"            angles - assigns min <float> and max <float> angle range where\n"
"                     selected angles should be in (e.g. angles\n"
"                     'chi1=90..160' or 'psi=60..70, phi=50..80');\n"
"            expand - expands the selection by one step that follows the\n"
"                     connections.\n"
"\n"
"        List operators:\n"
"            .. - range of integers (e.g. 2..5 is 2,3,4,5);\n"
"            ,  - list of integers or keywords (e.g. 2,3 or A,B).\n"
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
"        command (same as --target) describing the atoms that will be marked\n"
"        as selected (S) (default: target).\n"
"    --tags <tag>[,<tag>...]\n"
"        select PDBx tag that will be in the output\n"
"        (default: '_atom_site,\n"
"                   _rotag_rotamer_angle,\n"
"                   _rotag_dihedral_angle,\n"
"                   _rotag_rotamer_energy,\n"
"                   _rotag_pairwise_energy,\n"
"                   _rotag_energy,\n"
"                   _rotag_rmsd').\n"
"    -r, --related-data\n"
"        only related data records from other categories are shown when\n"
"        '_atom_site' records are selected.\n"
"    -p, --pdb\n"
"        input file is in PDB format.\n"
"    -k, --keep-ignored\n"
"        keep ignored atoms.\n"
"    -x, --rand-seed <int>\n"
"        set a seed for random (rand) selection.\n"
"    -v, --version\n"
"        print version" << std::endl;
                exit(1);
            case 'v':
                std::cout << version() << std::endl;
                exit(1);
        }
    }

    std::vector<std::string> tags_list;
    boost::replace_all(tags, " ", "");
    boost::split(tags_list, tags, boost::is_any_of(","));

    Parameters parameters(progname);

    // for (int index = optind; index < argc; index++) {
    //     AtomSite atom_site(argv[index], is_pdb);

    //     std::set<int64_t> target_atom_ids = selection_parser(
    //         parameters, atom_site, target_cmd, random_seed, 1
    //     );

    //     for (const int64_t& target_atom_id : target_atom_ids) {
    //         std::cout << target_atom_id << std::endl;
    //     }
    // }

    return 0;
}
