#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/leucine_014.cif

../programs/add_hydrogens ${pdbx_file}
