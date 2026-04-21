#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-selected-003.cif

rotag_scan --pairwise --limit 5 ${pdbx_file}
