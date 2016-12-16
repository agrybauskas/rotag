#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/tyrosine_017.cif

bond_length=1.592
bond_length_error=0.404

../programs/connect_atoms ${bond_length} ${bond_length_error} < ${cif_file} \
| sort
