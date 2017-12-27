#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/isoleucine_033.cif
target_resi="555"
angles="chi0 pi & chi1 pi & chi2 pi & chi3 pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
