#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_022.cif
target_residue="label_seq_id 18"

../programs/all_dihedral "${target_residue}" ${pdbx_file}
