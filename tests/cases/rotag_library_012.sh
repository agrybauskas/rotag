#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-006.cif

rotag_library -t 5 ${pdbx_file} 2>&1
