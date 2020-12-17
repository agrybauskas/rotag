#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-008.cif

rotag_scan \
    --parameters 'lj_k=0.0, c_k=0.0, h_k=1.0, t_k=0.0, cutoff_start=2.5, cutoff_end=5.0' \
    ${pdbx_file}
