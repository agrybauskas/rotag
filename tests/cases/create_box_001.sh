#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
atom_specifier="label_atom_id CA,C,CB,OG"
data_specifier="id,label_atom_id,label_comp_id,Cartn_x,Cartn_y,Cartn_z"

../programs/create_box "${atom_specifier}" "${data_specifier}" < ${cif_file}
