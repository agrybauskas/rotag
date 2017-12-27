#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/asparagine_020.cif

../programs/bond_type ${pdbx_file} | sort -k 1 -n
