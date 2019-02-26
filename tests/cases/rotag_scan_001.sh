#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_scan \
    --parameters 'lj_k=1.0, c_k=1.0, h_k=1.0, t_k=1.0, cutoff_start=2.5, cutoff_end=5.0' \
    --cutoff-atom Inf \
    ${pdbx_file}
