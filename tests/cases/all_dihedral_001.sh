#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_resi="label_seq_id 18"

../programs/all_dihedral "${target_resi}" ${cif_file}
