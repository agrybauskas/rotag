#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/leucine_034.cif
target_resi="546"
angles="chi0 0 & chi1 pi & chi2 pi & chi3 pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
