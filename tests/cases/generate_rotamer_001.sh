#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
target_resi="label_seq_id 219"
angle_range="chi0 0-0 & chi1 0-0".

../programs/generate_rotamer "${target_resi}" "${angle_range}" ${cif_file}
