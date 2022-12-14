#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

rotag_energy \
    --potential lennard_jones ${pdbx_file} \
    --parameters 'lj_k=1.0'
