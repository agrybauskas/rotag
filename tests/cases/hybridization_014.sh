#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/leucine_014.cif

../programs/hybridization ${pdbx_file}
