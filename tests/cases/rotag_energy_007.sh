#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_energy \
    --potential h_bond ${pdbx_file} \
    --parameters 'h_k=1.0'
