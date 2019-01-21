#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

rotag_energy \
    --potential composite ${pdbx_file} \
    --parameters 'lj_k=1.0, c_k=1.0, h_k=1.0, cutoff_start=2.5, cutoff_end=5.0'
