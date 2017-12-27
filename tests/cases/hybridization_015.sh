#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/methionine_015.cif

../programs/hybridization ${pdbx_file}
