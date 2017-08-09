#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif

../programs/obtain_atom_site ${pdbx_file}
