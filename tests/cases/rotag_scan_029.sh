#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-selected-002.cif

rotag_scan \
    -a 72.0 \
    ${pdbx_file}
