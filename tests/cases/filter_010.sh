#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
include_specifier="label_atom_id N,CA,C,CB,OG & label_comp_id SER"
exclude_specifier="type_symbol N,O & label_atom_id CA"
data_specifier="label_atom_id"

"$(dirname "$0")"/../scripts/filter "${include_specifier}" \
		                    "${exclude_specifier}" \
                		    "${data_specifier}"    \
		                    "${pdbx_dump_file}"
