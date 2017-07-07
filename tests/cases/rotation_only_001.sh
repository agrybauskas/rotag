#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"
angle_range="chi0 0-2*pi"
num_of_angles="chi0 16" # Number of angles witll be generated for each dihedral
                        # angle.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_angles}" \
			  ${cif_file}
