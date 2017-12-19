#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/threonine_019.cif

../programs/add_hydrogens ${pdbx_file}
