#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif
target_atom="label_atom_id NZ"
clash_specifier="is_pseudo_atom 1"
angle_range="chi0 0,2*pi & chi1 0,2*pi & chi2 0,2*pi & chi3 0,2*pi"
num_of_angles="chi0 2 & chi1 2 & chi2 2 & chi3 2" # Number of angles witll be
                                                  # generated for each dihedral
                                                  # angle.

../programs/radius_only "${target_atom}" "${clash_specifier}" "${angle_range}" \
			"${num_of_angles}" ${pdbx_file}
