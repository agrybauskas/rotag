#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/lysine_005.cif

../programs/add_hydrogen ${pdbx_file}
