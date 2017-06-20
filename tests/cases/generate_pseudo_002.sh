#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
target_atom="label_atom_id OD1"
angle_range="chi0 0-1*pi & chi1 0-1*pi"
num_of_points=12 # Number of pseudo-atoms that will be generated.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_points}" ${cif_file}
