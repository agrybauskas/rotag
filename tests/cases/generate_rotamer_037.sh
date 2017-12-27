#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/tyrosine_037.cif
target_resi="536"
angles="chi0 0.5*pi & chi1 0 & chi2 -0.5*pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
