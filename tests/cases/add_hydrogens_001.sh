#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/serine_001.cif

../programs/add_hydrogens ${pdbx_file}
