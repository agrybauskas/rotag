#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id CA,N,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"
number_of_atoms=17 # Number of pseudo-atoms that will be generated.

../programs/bond_torsion "${atom_specifier}" "${data_specifier}" \
			 "${number_of_atoms}" ${pdbx_file}
