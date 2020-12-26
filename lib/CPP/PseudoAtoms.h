#ifndef _PSEUDOATOMS_H_
#define _PSEUDOATOMS_H_

#include "PseudoAtoms.h"

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <map>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

void generate_library(SV* hash_ref);

#endif
