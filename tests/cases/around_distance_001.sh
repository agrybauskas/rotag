#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
include_specifier="label_atom_id CA"
around_distance=2.41

"$(dirname "$0")"/../scripts/around_distance "${include_specifier}" \
                                             ${around_distance} \
		                             "${pdbx_dump_file}"
