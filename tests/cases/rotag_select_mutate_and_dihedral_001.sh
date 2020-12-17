#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

rotag_select ${pdbx_file} \
    | rotag_mutate -m '1:SER,chi1=0.0' \
    | rotag_dihedral -S
