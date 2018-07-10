#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-008.cif

rotag_library -i ${pdbx_file} --parameters 'lj_epsilon 0; c_k 0' --top-rank 1
