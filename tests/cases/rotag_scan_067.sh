#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-004.cif

rotag_scan --use-library -H ${pdbx_file}
