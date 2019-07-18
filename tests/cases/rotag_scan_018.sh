#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_scan \
    --parameters 'lj_k=0.0, c_k=0.0, h_k=1.0, t_k=0.0, cutoff_start=2.5, cutoff_end=5.0' \
    --categories '_[local]_rotamer_angle' \
    ${pdbx_file} 2>&1
