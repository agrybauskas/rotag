#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif
target_resi="219"
angles="chi0 0 & chi1 0"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
