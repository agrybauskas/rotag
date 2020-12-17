#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/sidechain-stream-001.cif

rotag_rmsd -d -c '1,2' --tags '_[local]_rmsd' -F csv ${pdbx_file}
