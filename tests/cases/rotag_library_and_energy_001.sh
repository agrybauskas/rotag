#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -c 'Inf' -t 1 ${pdbx_file} \
    | rotag_energy -r
