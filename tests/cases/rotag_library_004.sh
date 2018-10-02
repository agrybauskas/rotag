#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-004.cif

rotag_library --parameters 'lj_epsilon = 0, h_epsilon = 0' --top-rank 1 ${pdbx_file}
