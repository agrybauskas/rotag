#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamine_021.cif
target_resi="5"
angles="chi0 pi & chi1 0 & chi2 pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
