#!/bin/bash

pdbx_file_1=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_rmsd -c '2,3' ${pdbx_file_1} ${pdbx_file_2} 2>&1
