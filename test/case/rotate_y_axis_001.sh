#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
angle="1.5708" # In radians.

../programs/rotate_y_axis "${atom_specifier}" "${data_specifier}" \
			  "${angle}" < ${cif_file}
