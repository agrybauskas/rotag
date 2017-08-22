#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
angle_min_max=1.0 # +- value that will change in bond angle (radians).
num_of_atoms=27 # Number of pseudo-atoms that will be generated.

../programs/angle_bending "${atom_specifier}" "${data_specifier}" \
			  "${angle_min_max}" "${num_of_atoms}" \
			  ${pdbx_file}
