#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_rmsd -c '1,2,3' ${pdbx_file} 2>&1
