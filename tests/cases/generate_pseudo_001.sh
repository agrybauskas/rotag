#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"
angle_range="chi0 1*pi-1*pi"
num_of_angles="chi0 1" # Number of angles will be generated for each dihedral
                       # angle.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_angles}" ${cif_file}
