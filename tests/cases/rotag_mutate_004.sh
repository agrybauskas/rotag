#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-selected-002.cif
moiety_file=$(dirname "$0")/../inputs/moieties/arginine-moiety.cif

rotag_mutate \
    -m 'ARG' \
    -M ${moiety_file} \
    -a 'chi1=90.0, chi2=0.0, chi3=90.0, chi4=0.0' \
    ${pdbx_file}
