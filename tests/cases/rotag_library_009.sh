#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/lysine-library-001.cif

rotag_library -t 100 ${pdbx_file}
