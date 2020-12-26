#include "PseudoAtoms.h"

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <map>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
/* --------------------- Contructors and/or destructors ---------------------- */

void generate_library(SV* hash_ref)
{
    HV * hash = (HV*) SvRV( hash_ref );
    SV** value = hv_fetch(hash, "a", 1, 0);
    std::cout << SvPV_nolen(*value) << std::endl;
}

/* -------------------------------- Methods ---------------------------------- */

/* -------------------------- Setters and Getters ---------------------------- */
