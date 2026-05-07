#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-selected-004.cif

rotag_scan -H --pairwise --limit 5 -x 23 -X 5 ${pdbx_file}
