#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_rmsd -S -d --tags '_[local]_rmsd,_[local]_atom_site' -F csv -c '1,2' ${pdbx_file} 2>&1
