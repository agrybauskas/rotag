#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/isoleucine_013.cif

../programs/connect_atoms ${cif_file} | sort -k 1 -n
