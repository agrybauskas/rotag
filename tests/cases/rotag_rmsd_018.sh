#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamine-002.cif

rotag_rmsd -c '2,4' ${pdbx_file} 2>&1
