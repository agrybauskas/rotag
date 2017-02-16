#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

angle_min_max=10 # +- value that will change in bond length.
num_of_points=100 # Number of points that will be generated.

../programs/change_bond_angle "${atom_specifier}" "${data_specifier}" \
			      "${angle_min_max}" "${num_of_points}" < ${cif_file}
