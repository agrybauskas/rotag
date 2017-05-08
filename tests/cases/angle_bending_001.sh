#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

angle_min_max=1.0 # +- value that will change in bond length.
num_of_points=27 # Number of points that will be generated.

../programs/angle_bending "${atom_specifier}" "${data_specifier}" \
			  "${angle_min_max}" "${num_of_points}" \
			  ${cif_file}
