#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_rmsd -d -c '1,2' -F csv ${pdbx_file} 2>&1
