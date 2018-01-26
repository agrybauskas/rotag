#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/methionine_015.cif
target_resi="523"
angles="chi0 pi & chi1 -1*pi & chi2 0.5*pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
