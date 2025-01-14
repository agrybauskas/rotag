#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-with-h2o-and-mg-model-001.cif

rotag_rmsd -b -S -c '1,2' ${pdbx_file}
