#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-002.cif

rotag_rmsd -S -d --no-strict --related-data -c  '1,2' ${pdbx_file}
