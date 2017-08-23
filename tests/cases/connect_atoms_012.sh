#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/valine_012.cif

../programs/connect_atoms ${pdbx_file} | sort -k 1 -n
