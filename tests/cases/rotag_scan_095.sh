#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/glutamic-acid-with-h2o-and-mg-003.cif

rotag_scan -H --rand-seed 23 --rand-step 10 --top-rank 10 ${pdbx_file}
