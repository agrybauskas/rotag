#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-006.cif

rotag_rmsd -S -c  '1,2' ${pdbx_file}
