#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-selected-002.cif

rotag_rmsd -c '1,2' ${pdbx_file} 2>&1
