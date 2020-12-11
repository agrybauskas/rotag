#!/bin/bash

pdbx_file_1=$(dirname "$0")/../inputs/models/aspartic-acid-model-010.cif
pdbx_file_2=$(dirname "$0")/../inputs/models/aspartic-acid-model-011.cif

rotag_rmsd -b -S -c  '1,2' ${pdbx_file_1} ${pdbx_file_2}
