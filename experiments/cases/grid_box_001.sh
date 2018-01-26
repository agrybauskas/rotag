#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id N,CA,CB,OG,C,O"

../programs/grid_box "${atom_specifier}" ${pdbx_file}
