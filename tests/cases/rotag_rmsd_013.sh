#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_rmsd -S ${pdbx_file} 2>&1
