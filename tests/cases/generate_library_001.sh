#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif
target_resi="18"

../programs/generate_library "${target_resi}" ${pdbx_file}
