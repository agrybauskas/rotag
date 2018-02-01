#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamine_021.cif
start_atom="label_atom_id CA"
next_atom="label_atom_id CB"

../programs/rotatable_bonds "${start_atom}" "${next_atom}" ${pdbx_file}
