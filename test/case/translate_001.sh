#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
transl_coord="1.000,1.000,1.000"

../programs/translate "${atom_specifier}" "${data_specifier}" \
		      "${transl_coord}" < ${cif_file}
