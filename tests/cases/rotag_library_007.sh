#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-007.cif

rotag_library \
    --parameters 'lj_epsilon=0.0,c_k=0.0,h_epsilon=1.0,r_sigma=2.0,cutoff_start=2.5,cutoff_end=5.0' \
    --top-rank 1 \
    ${pdbx_file}
