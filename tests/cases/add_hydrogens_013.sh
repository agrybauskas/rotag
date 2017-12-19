#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/isoleucine_013.cif

../programs/add_hydrogens ${pdbx_file}
