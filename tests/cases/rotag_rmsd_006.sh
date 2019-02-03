#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/sidechain-stream-001.cif

rotag_rmsd -d -c '1,2' -F csv ${pdbx_file}
