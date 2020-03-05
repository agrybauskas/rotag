#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/tyrosine-001.cif

rotag_scan -a 360.00 -p 'cutoff_atom=Inf' ${pdbx_file} \
    | rotag_add -S
