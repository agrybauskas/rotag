#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -c 'Inf' -C 'Inf' -t 5 ${pdbx_file} \
    | rotag_add -r \
    | rotag_select \
    | rotag_dihedral
