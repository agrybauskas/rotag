#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi"
target_atom="label_atom_id OD1"
num_of_angles="chi0 4 & chi1 4"  # Number of angles witll be generated for each
                                 # dihedral angle.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_angles}" \
			  ${cif_file}
