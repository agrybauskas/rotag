#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamic_acid_007.cif
target_residue="label_seq_id 14"

../programs/all_dihedral "${target_residue}" ${pdbx_file}
