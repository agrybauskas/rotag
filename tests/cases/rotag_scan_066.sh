#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-003.cif

rotag_scan --use-library -H ${pdbx_file}
