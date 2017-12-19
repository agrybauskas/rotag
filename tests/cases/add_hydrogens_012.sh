#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/valine_012.cif

../programs/add_hydrogens ${pdbx_file}
