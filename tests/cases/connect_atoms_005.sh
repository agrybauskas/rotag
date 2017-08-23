#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif

../programs/connect_atoms ${pdbx_file} | sort -k 1 -n
