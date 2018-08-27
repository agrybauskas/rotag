#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -c 'Inf' -C 'Inf' ${pdbx_file}
