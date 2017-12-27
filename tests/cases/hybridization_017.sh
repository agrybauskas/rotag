#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/tyrosine_017.cif

../programs/hybridization ${pdbx_file}
