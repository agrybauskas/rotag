#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_021.cif
atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

../programs/check_ref_frame "${atom_specifier}" "${data_specifier}" ${pdbx_file}
