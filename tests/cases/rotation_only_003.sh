#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamic_acid_007.cif
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi"
atom_specifier="label_atom_id OE1"
num_of_angles="chi0 3 & chi1 3 & chi2 3"  # Number of angles witll be generated
                                          # for each dihedral angle.

../programs/rotation_only "${atom_specifier}" "${angle_range}" \
			  "${num_of_angles}" ${pdbx_file}
