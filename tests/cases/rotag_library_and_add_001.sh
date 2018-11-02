#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library ${pdbx_file} -t 1 \
    | rotag_add -r
