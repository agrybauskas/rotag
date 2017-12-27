#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/cysteine_008.cif

../programs/hybridization ${pdbx_file}
