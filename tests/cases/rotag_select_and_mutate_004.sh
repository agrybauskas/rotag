#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/5xjp.cif

rotag_select -t '1:resid 20,21' ${pdbx_file}  | \
    rotag_mutate -m '1:PHE,chi1=0.0,chi2=0.0'
