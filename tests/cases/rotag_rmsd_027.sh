#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-007.cif

rotag_rmsd -S -b -c  '1,2' ${pdbx_file}
