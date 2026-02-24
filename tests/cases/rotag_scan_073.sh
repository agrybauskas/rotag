#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-003.cif

rotag_scan --use-library -H -M 0.1 ${pdbx_file}
