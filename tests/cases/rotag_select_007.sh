#!/bin/bash

pdbx_file1=$(dirname "$0")/../inputs/5svd.cif
pdbx_file2=$(dirname "$0")/../inputs/5aox.cif

rotag_select -t 'atomid 4205' -s 'atomid 4206' \
             -t 'atomid 57'   -s 'atomid 58' \
             ${pdbx_file1}    ${pdbx_file2}
