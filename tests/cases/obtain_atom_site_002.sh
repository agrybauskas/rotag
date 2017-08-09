#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/5svd_002.cif

../programs/obtain_atom_site ${pdbx_file}
