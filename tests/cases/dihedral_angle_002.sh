#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_021.cif
atom_specifier="label_atom_id C,O,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

../programs/dihedral_angle "${atom_specifier}" "${data_specifier}" ${pdbx_file}
