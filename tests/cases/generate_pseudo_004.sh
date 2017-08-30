#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif
target_atom="label_atom_id NZ"
angle_range="chi0 1.5*pi,1.5*pi & chi1 1.5*pi,1.5*pi & chi2 1.5*pi,1.5*pi\
             & chi3 1.5*pi,1.5*pi"
num_of_angles="chi0 1 & chi1 1 & chi2 1 & chi3 1" # Number of angles will be
                                                  # generated for each dihedral
                                                  # angle.

../programs/generate_pseudo "${target_atom}" "${angle_range}" \
			    "${num_of_angles}" ${pdbx_file}
