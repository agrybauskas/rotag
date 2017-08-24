#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif
target_residue="label_seq_id 219"

../programs/all_dihedral "${target_residue}" ${pdbx_file}
