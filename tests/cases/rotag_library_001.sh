#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-001.cif

rotag_library -t 1 ${pdbx_file}
