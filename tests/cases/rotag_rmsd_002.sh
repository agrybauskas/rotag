#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamine-002.cif

rotag_rmsd -c '1,2' -F csv ${pdbx_file}
