#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif
target_residue="label_seq_id 572"

../programs/all_dihedral "${target_residue}" ${pdbx_file}
