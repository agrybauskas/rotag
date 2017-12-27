#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/valine_032.cif
target_resi="100"
angles="chi0 0 & chi1 0 & chi2 0"

../programs/generate_rotamer "${target_resi}" "${angles}" ${pdbx_file}
