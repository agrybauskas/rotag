#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/5xjp.cif

rotag_select -t '1:resid 32,33,34 && chain A' -s 'byres(target around 5) && ! hetatoms' -k ${pdbx_file}  | \
    rotag_mutate -m '1:TRP'
