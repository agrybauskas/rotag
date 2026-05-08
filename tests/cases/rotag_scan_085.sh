#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-selected-006.cif

rotag_scan -H -x 23 -X 2 --pairwise --limit 2 ${pdbx_file}
