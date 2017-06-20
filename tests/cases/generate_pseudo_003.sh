#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif
target_atom="label_atom_id OE1"
angle_range="chi0 0-1*pi & chi1 0-1*pi & chi2 0-1*pi"
num_of_points=16 # Number of pseudo-atoms that will be generated.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_points}" ${cif_file}
