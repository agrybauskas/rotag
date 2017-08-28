#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif
angle_range="chi0 -pi,pi & chi1 -pi,pi & chi2 -pi,pi & chi3 -pi,pi"
atom_specifier="label_atom_id NZ"
num_of_angles="chi0 2 & chi1 2 & chi2 2 & chi3 2" # Number of angles witll be
                                                  # generated for each dihedral
                                                  # angle.

../programs/rotation_only "${atom_specifier}" "${angle_range}" \
			  "${num_of_angles}" ${pdbx_file}
