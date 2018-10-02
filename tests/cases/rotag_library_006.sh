#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-006.cif

rotag_library --parameters 'lj_epsilon = 0, c_k = 0' --top-rank 1 ${pdbx_file}
