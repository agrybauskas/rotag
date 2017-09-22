#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"
clash_specifier="is_pseudo_atom 1"
angle_range="chi0 0,2*pi"
num_of_angles="chi0 16" # Number of pseudo-atoms that will be generated.

../programs/radius_only "${target_atom}" "${clash_specifier}" "${angle_range}" \
			"${num_of_angles}" ${pdbx_file}
