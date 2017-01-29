#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,C,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

rotational_angle="2*\$pi"

../programs/rotate_bond "${atom_specifier}" "${data_specifier}" "${rotational_angle}" < ${cif_file}
