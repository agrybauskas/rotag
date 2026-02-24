#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-002.cif

rotag_scan --use-library -H -M 0.05 ${pdbx_file}
