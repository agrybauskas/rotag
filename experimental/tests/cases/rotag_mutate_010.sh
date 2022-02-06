#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/cysteine-001.cif

rotag_mutate -m '1:SER,,chi1=3.14' ${pdbx_file} 2>&1
