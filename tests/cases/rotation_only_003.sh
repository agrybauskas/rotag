#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi"
target_atom="label_atom_id OE1"
num_of_angles="chi0 3 & chi1 3 & chi2 3"  # Number of angles witll be generated
                                          # for each dihedral angle.

../programs/rotation_only "${target_atom}" "${angle_range}" "${num_of_angles}" \
			  ${cif_file}
