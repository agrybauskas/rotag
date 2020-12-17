#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/phenylalanine-library-001.cif

rotag_library -F csv 2>&1 ${pdbx_file}
