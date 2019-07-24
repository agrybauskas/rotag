#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

cat ${pdbx_file} \
    | rotag_energy \
          --potential lennard_jones  \
          --pairwise \
          --parameters 'lj_k=1.0'
