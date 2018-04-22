#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

matrices_file=$(dirname "$0")/../inputs/matrices/matrices-003.dat
variable_values="x=6,y=3,z=2"

$(dirname "$0")/../scripts/mult_matrix_product "${variable_values}" \
	                                       "${matrices_file}" 2>&1 \
| sed 's/line\s*[0-9]*.$/line <row>./g'
