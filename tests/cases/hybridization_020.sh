#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/asparagine_020.cif

../programs/hybridization ${pdbx_file}
