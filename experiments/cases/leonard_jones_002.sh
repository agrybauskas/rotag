#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif
target_atom="label_atom_id OD1"
clash_specifier="is_pseudo_atom 1"
cutoff=0.04
angle_range="chi0 0,2*pi & chi1 0,2*pi"
num_of_angles="chi0 4 & chi1 4" # Number of pseudo-atoms that will be generated.

../programs/leonard_jones "${target_atom}" "${clash_specifier}" "${cutoff}" \
			  "${angle_range}" "${num_of_angles}" ${pdbx_file}
