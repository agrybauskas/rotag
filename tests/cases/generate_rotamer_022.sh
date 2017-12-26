#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_022.cif
target_resi="18"
angles="chi0 -0.5*pi & chi1 0"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
