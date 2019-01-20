#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.cif

rotag_library -a 1.00 -c Inf ${pdbx_file} \
    | rotag_add -S
