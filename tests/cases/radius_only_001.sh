#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"
angle_range="chi0 0-2*pi"
num_of_points=16 # Number of pseudo-atoms that will be generated.

../programs/radius_only "${target_atom}" "${angle_range}" "${num_of_points}" \
			  ${cif_file}
