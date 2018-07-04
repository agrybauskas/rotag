#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-004.cif

rotag_library -i ${pdbx_file} --parameters 'lj_epsilon 0; h_epsilon 0' --top-rank 1
