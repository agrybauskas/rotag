#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glycine_009.cif

../programs/hybridization ${pdbx_file}
