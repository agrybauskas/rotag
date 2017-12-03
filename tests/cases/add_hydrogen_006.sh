#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif

../programs/add_hydrogen ${pdbx_file}
