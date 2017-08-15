#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
angle="1.5708" # In radians.

../programs/x_axis_rotation "${atom_specifier}" "${data_specifier}" \
			    "${angle}" ${pdbx_file}
