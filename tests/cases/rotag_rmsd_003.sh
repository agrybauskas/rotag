#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamine-002.cif

rotag_rmsd -c '1,2; 3,4' --tags '_[local]_rmsd' -F csv ${pdbx_file}
