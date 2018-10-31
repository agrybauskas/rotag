#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_energy \
    --potential h_bond ${pdbx_file} \
    --parameters 'h_epsilon=1.0, r_sigma=2.0, cutoff_atom=0.5, cutoff_residue=1.0'
