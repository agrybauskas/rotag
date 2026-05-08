#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-selected-005.cif

rotag_scan --pairwise --limit 2 ${pdbx_file}
