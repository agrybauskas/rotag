#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/isoleucine_013.cif

../programs/hybridization ${pdbx_file}
