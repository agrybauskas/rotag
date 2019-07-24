#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/phenylalanine-library-001.cif

cat ${pdbx_file} | rotag_library -t 1
