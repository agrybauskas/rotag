#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif
target_atom="label_atom_id OE1"
angle_range="chi0 0.5*pi-0.5*pi & chi1 0.5*pi-0.5*pi & chi2 0.5*pi-0.5*pi"
num_of_angles="chi0 1 & chi1 1 & chi2 1" # Number of angles witll be generated
                                         # for each dihedral angle.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_angles}" ${cif_file}
