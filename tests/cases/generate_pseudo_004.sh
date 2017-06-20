#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/lysine_005.cif
target_atom="label_atom_id NZ"
angle_range="chi0 0-1*pi & chi1 0-1*pi & chi2 0-1*pi & chi3 0-1*pi"
num_of_points=8 # Number of pseudo-atoms that will be generated.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_points}" ${cif_file}
