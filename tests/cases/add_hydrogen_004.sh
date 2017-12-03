#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/phenylalanine_016.cif

../programs/add_hydrogen ${pdbx_file}
