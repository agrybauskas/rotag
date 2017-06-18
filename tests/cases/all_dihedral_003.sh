#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif
target_resi="label_seq_id 14"

../programs/all_dihedral "${target_resi}" ${cif_file}
