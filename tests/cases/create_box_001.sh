#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
atom_specifier="label_atom_id CA,C,CB,OG"

../programs/create_box "${atom_specifier}" ${pdbx_file}
