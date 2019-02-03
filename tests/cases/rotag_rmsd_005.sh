#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamine-002.cif

rotag_rmsd -d -c '1,2; 3,4' ${pdbx_file}
