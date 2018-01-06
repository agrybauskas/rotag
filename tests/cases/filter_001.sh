#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
include_specifier="label_atom_id CA,C,CB,OG"
exclude_specifier="label_atom_id CA"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

../programs/filter "${include_specifier}" \
		   "${exclude_specifier}" \
		   "${data_specifier}"    \
		   ${pdbx_file}
