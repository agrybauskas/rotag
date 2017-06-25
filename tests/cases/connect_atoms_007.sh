#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/glutamic_acid_007.cif

../programs/connect_atoms ${cif_file} | sort -k 1 -n
