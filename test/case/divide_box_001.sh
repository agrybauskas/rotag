#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif
atom_specifier="label_atom_id CA,C,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

../programs/divide_box "${atom_specifier}" "${data_specifier}" < ${cif_file}
