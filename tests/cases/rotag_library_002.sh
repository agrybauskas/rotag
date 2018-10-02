#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_library --parameters 'c_k = 0, h_epsilon = 0' --top-rank 1 ${pdbx_file}
