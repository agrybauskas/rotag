#include "PseudoAtoms.h"

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <map>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

void generate_library(SV* options_ptr)
{
    HV* options = (HV*) SvRV(options_ptr);

    HV* atom_site = (HV*) SvRV((SV*) hv_fetch(options, "atom_site", 9, 0));
    HV* residue_unique_keys =
        (HV*) SvRV((SV*) hv_fetch(options, "residue_unique_keys", 19, 0));
    HV* include_interactions =
        (HV*) SvRV((SV*) hv_fetch(options, "include_interactions", 20, 0));
    HV* angles = (HV*) SvRV((SV*) hv_fetch(options, "angles", 6, 0));
    // NOTE: True and false comes from key existance.
    bool rmsd =
        hv_exists(options,"rmsd",4) && SvUV(*hv_fetch(options, "rmsd",4,0)) != 0;
    HV* interactions = (HV*) SvRV((SV*) hv_fetch(options, "interactions", 12,0));
    AV* interaction_atom_names =
        (AV*) SvRV((SV*) hv_fetch(options, "interaction_atom_names", 22,0));
    double edge_length_interaction =
        SvNV(*hv_fetch(options, "edge_length_interaction", 23,0));
    // TODO: For now, threads are not included, but should be in the future.
    HV* suboptions = (HV*) SvRV((SV*) hv_fetch(options, "options", 7, 0));

    // HACK: the variables should be changeable in the future.
    std::string conf_model =  "rotation_only";
    int threads = 1;
    // include_interactions

    // potential_functions

    // rotamer_library

    // atom_site_groups

    // HV * options = (HV*) SvRV( options_ptr );
    // SV** value = hv_fetch(hash, "a", 1, 0);
    // std::cout << SvPV_nolen(*value) << std::endl;
}
