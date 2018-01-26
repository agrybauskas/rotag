#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/methionine_035.cif
target_resi="523"
angles="chi0 pi & chi1 -1*pi & chi2 0.5*pi & chi3 0.5*pi & chi4 0.5*pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
