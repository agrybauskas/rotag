#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/lysine_005.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi & chi3 0-2*pi"
target_atom="label_atom_id NZ"
num_of_points=8 # Number of pseudo-atoms that will be generated.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_points}" \
			  ${cif_file}
