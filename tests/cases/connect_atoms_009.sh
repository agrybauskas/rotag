#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glycine_009.cif

bond_length=1.592
bond_length_error=0.404

../programs/connect_atoms ${bond_length} ${bond_length_error} ${cif_file} \
| sort -k 1 -n
