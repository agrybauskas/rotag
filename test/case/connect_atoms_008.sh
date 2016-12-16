#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/cysteine_008.cif

bond_length=1.360
bond_length_error=0.174

../programs/connect_atoms ${bond_length} ${bond_length_error} < ${cif_file} \
| sort
