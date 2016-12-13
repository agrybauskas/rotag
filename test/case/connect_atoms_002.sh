#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/5svd_002.cif
atom_specifier="label_atom_id CA,C,CB,OG"
data_specifier="Cartn_x,Cartn_y,Cartn_z"

../programs/connect_atoms "${atom_specifier}" "${data_specifier}" < ${cif_file}
