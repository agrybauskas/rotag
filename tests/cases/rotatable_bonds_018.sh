#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/tryptophan_018.cif
start_atom="label_atom_id CA"

../programs/rotatable_bonds "${start_atom}" ${pdbx_file}
