#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
include_specifier=
exclude_specifier="label_atom_id CA,C,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

"$(dirname "$0")"/../scripts/filter "${include_specifier}" \
		                    "${exclude_specifier}" \
                		    "${data_specifier}"    \
		                    "${pdbx_dump_file}"
