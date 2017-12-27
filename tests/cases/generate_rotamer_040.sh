#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/asparagine_040.cif
target_resi="135"
angles="chi0 pi & chi1 0"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
