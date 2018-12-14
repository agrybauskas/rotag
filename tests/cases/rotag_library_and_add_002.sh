#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_library ${pdbx_file} -c Inf \
    | rotag_add -r
