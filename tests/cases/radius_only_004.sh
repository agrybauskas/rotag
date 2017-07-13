#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/lysine_005.cif
target_atom="label_atom_id NZ"
clash_specifier="is_pseudo_atom 1"
angle_range="chi0 0-2*pi & chi1 0-2*pi & chi2 0-2*pi & chi3 0-2*pi"
num_of_angles="chi0 4 & chi1 4 & chi2 4 & chi3 4" # Number of angles witll be
                                                  # generated for each dihedral
                                                  # angle.

../programs/radius_only "${target_atom}" "${clash_specifier}" "${angle_range}" \
			"${num_of_angles}" ${cif_file}
