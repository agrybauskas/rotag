#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/histidine_004.cif

../programs/hybridization ${pdbx_file}
