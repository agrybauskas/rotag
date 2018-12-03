#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-001.cif

rotag_library \
    -H \
    -c 40.0 \
    -C 40.0 \
    --parameters 'lj_epsilon=1.0,c_k=1.0,h_epsilon=1.0,r_sigma=2.0,cutoff_start=2.5,cutoff_end=5.0' \
    --top-rank 1 \
    ${pdbx_file}
