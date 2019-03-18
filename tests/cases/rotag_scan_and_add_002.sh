#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_scan ${pdbx_file} -p 'cutoff_atom=Inf' \
    | rotag_add -S
