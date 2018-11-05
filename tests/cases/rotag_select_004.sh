#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/2cj3.cif

rotag_select -x 22 -t 'byres((atomname CA & resname TRP) rand 2)' ${pdbx_file} 2>&1 \
    | sed 's/line\s*[0-9]*/line <row>/g'
