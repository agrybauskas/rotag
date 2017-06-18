#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/lysine_005.cif
target_resi="label_seq_id 572"

../programs/all_dihedral "${target_resi}" ${cif_file}
