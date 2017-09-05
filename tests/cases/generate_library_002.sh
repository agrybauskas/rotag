#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif
target_resi="219"
small_angle="0.1*pi"

../programs/generate_library "${target_resi}" "${small_angle}" ${pdbx_file}
