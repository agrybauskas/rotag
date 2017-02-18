#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
ref_frame="global"

../programs/switch_ref_frame "${atom_specifier}" "${data_specifier}" \
			     "${ref_frame}" < ${cif_file}
