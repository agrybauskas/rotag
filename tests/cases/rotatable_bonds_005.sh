#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif
start_atom="label_atom_id CA"
next_atom="label_atom_id CB"

../programs/rotatable_bonds "${start_atom}" "${next_atom}" ${pdbx_file}
