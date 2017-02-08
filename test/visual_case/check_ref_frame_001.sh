#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

num_of_points=17 # Number of points that will be generated

../programs/check_ref_frame "${atom_specifier}" "${data_specifier}" \
			    "${num_of_points}" < ${cif_file}
