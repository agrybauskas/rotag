#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/arginine_003.cif

../programs/add_hydrogens ${pdbx_file}
