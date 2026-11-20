#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/glutamic-acid-with-h2o-and-mg-005.cif

rotag_scan -H --rand-seed 23 --rand-step 5 --top-rank 5 ${pdbx_file}
