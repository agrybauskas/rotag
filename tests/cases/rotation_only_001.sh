#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id OG"
angle_range="chi0 0*pi,2*pi"
num_of_angles="chi0 16" # Number of angles witll be generated for each dihedral
                        # angle.

../programs/rotation_only "${atom_specifier}" "${angle_range}" \
			  "${num_of_angles}" ${pdbx_file}
