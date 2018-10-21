#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library \
    --parameters 'lj_epsilon=1.0,c_k=1.0,h_epsilon=1.0,r_sigma=2.0,cutoff_start=2.5,cutoff_end=5.0' \
    --cutoff-atom 0.4 \
    --cutoff-residue 1.0 \
    ${pdbx_file}
