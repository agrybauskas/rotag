#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

rotag_energy \
    --potential lennard_jones ${pdbx_file} \
    --parameters 'lj_epsilon=1.0, cutoff_atom=0.5, cutoff_residue=1.0'
