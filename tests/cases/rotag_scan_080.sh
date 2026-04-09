#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-005.cif

rotag_scan --use-library -x 23 -X 10 ${pdbx_file}
