#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-amino-acids-with-connections-library-003.cif

rotag_scan --use-library -H --rand-seed 23 --rand-step 5 --top-rank 5 ${pdbx_file}
