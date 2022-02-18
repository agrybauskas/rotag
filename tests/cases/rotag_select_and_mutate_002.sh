#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/5xjp.cif

rotag_select -t '1:resid 20,21 && chain A' -s 'byres(target around 5) && ! hetatoms' ${pdbx_file}  | \
    rotag_mutate -m '1:TYR'
