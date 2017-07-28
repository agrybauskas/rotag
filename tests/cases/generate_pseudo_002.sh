#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
target_atom="label_atom_id OD1"
angle_range="chi0 0-0 & chi1 0-0"
num_of_angles="chi0 1 & chi1 1" # Number of angles will be generated for each dihedral
                 # angle.
../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_angles}" ${cif_file}
