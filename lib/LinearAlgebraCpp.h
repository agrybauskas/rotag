#ifndef _LINEARALGEBRACPP_H_
#define _LINEARALGEBRACPP_H_

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

void create_ref_frame( SV*, SV*, SV* );

double calculate_vector_length( double* );

#endif
