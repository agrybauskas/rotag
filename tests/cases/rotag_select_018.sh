#!/bin/bash

pdbx_file_1=$(dirname "$0")/../inputs/amino-acids/asparagine-001.cif
pdbx_file_2=$(dirname "$0")/../inputs/amino-acids/lysine-001.cif

rotag_select \
    -t 'atomname C' -t 'atomname N' -s 'atomname CA'  \
    ${pdbx_file_1} ${pdbx_file_2}
