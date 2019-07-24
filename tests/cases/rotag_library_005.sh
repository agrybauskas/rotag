#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/phenylalanine-library-001.cif

rotag_library -T 0.001 ${pdbx_file}
