#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/asparagine_020.cif

../programs/add_hydrogens ${pdbx_file}
