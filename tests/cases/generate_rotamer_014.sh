#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/leucine_014.cif
target_resi="546"
angles="chi0 0 & chi1 0.5*pi"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
