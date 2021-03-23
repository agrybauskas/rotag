#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_scan -a '-180.0..36.0..180.0' -p 'cutoff_atom=Inf' ${pdbx_file} \
    | rotag_add -S \
    | rotag_select \
    | rotag_dihedral -S -r
