#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamic_acid_007.cif
target_resi="14"
angles="chi0 0 & chi1 pi & chi2 0"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
