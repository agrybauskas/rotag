#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/lysine_005.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi & chi3 0-2*pi"
target_atom="label_atom_id NZ"
num_of_angles="chi0 2 & chi1 2 & chi2 2 & chi3 2" # Number of angles witll be
                                                  # generated for each dihedral
                                                  # angle.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_angles}" \
			  ${cif_file}
