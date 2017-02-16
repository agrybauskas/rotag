#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

bond_min_max=2.0 # +- value that will change in bond length.
num_of_points=10 # Number of points that will be generated.

../programs/change_bond_len "${atom_specifier}" "${data_specifier}" \
			    "${bond_min_max}" "${num_of_points}" < ${cif_file}
