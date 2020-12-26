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
    double edge_length_interation;
    // interaction_atom_names;
    // cutoff_atom;

    // conf_model
    // threads
    // include_interactions

    // potential_functions

    // rotamer_library

    // atom_site_groups

    // HV * options = (HV*) SvRV( options_ptr );
    // SV** value = hv_fetch(hash, "a", 1, 0);
    // std::cout << SvPV_nolen(*value) << std::endl;
}
