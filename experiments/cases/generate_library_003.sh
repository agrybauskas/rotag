#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamic_acid_007.cif
target_resi="14"
small_angle="0.1*pi"

../programs/generate_library "${target_resi}" "${small_angle}" ${pdbx_file}
