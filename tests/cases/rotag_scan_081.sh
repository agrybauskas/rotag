#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-005.cif

rotag_scan --use-library -H -x 23 -X 10 --limit 10 ${pdbx_file}
