#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamic_acid_007.cif

../programs/add_hydrogens ${pdbx_file}
