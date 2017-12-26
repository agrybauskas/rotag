#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/isoleucine_013.cif
target_resi="555"
angles="chi0 pi & chi1 0.5*pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
