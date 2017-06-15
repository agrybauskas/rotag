#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
target_resi="label_seq_id 219"

../programs/all_dihedral "${target_resi}" ${cif_file}
