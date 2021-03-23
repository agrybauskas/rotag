#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

rotag_scan \
    --rand-seed 23 \
    --rand-count 1 \
    --parameters 'cutoff_atom=Inf, lj_k=0.0, c_k=1.0, h_k=0.0, t_k=0.0, cutoff_start=2.5, cutoff_end=5.0' \
    ${pdbx_file}
