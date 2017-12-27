#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/alanine_011.cif

../programs/hybridization ${pdbx_file}
