#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.cif

rotag_scan -a 1.00 -p 'cutoff_atom=Inf' ${pdbx_file} \
    | rotag_add -S
