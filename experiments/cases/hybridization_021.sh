#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/glutamine_021.cif

../programs/hybridization ${pdbx_file}
