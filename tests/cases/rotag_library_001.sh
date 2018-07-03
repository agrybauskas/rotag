#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -i ${pdbx_file} --potential leonard_jones --parameters 'c_k 0; h_epsilon 0' --cutoff-atom 10 --cutoff-residue 0
