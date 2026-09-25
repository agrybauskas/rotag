#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-009.cif

rotag_scan --verbose --verbosity-level 3 -H --rand-seed 23 --rand-step 5 --top-rank 5 ${pdbx_file} 2>&1
