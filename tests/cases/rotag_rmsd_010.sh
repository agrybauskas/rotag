#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_rmsd -S -b -d --related-data -c '1,2' ${pdbx_file}
