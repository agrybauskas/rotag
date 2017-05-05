#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi"
target_atom="label_atom_id OE1"
num_of_points=16 # Number of pseudo-atoms that will be generated.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_points}" \
			  < ${cif_file}
