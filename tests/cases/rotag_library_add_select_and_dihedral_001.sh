#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -c 'Inf' --top-rank 5 ${pdbx_file} \
    | rotag_add -S \
    | rotag_select \
    | rotag_dihedral -S -r
