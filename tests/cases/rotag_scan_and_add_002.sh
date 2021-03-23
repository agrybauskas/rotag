#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_scan -a '-180.0..36.0..180.0' ${pdbx_file} -p 'cutoff_atom=Inf' \
    | rotag_add -S
