#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
target_resi="18"
small_angle="0.1*pi"

../programs/generate_library "${target_resi}" "${small_angle}" ${pdbx_file}
